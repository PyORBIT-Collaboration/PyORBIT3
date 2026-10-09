import math

import pytest

from orbit.core.spacecharge import (
    Grid2D,
    PoissonSolverFFT2D,
    SpaceChargeCalc2p5D,
    SpaceChargeCalc2p5Drb,
    SpaceChargeCalcSliceBySlice2D,
)


def _antiderivative(x, y):
    if x == 0.0 or y == 0.0:
        return 0.0
    return (
        -x * y * math.log(math.hypot(x, y))
        + 1.5 * x * y
        - x * x * math.atan(y / x) / 2.0
        - y * y * math.atan(x / y) / 2.0
    )


def _cell_average(dx, dy, x, y):
    def f(ix, iy):
        return _antiderivative((ix - 0.5) * dx, (iy - 0.5) * dy)

    return (f(x + 1, y + 1) - f(x, y + 1) - f(x + 1, y) + f(x, y)) / (dx * dy)


def _solve(use_integrated=False):
    shape = (9, 8)
    limits = (-1.0, 1.0, -2.0, 2.0)
    solver = PoissonSolverFFT2D(*shape, *limits)
    solver.setUseIntegratedGreenFunction(use_integrated)
    rho = Grid2D(*shape, *limits)
    phi = Grid2D(*shape, *limits)
    source = (4, 4)
    rho.setValue(1.0, *source)
    solver.findPotential(rho, phi)
    return solver, phi, source


def test_integrated_kernel_matches_direct_cell_integral_on_rectangular_grid():
    _, phi, source = _solve(use_integrated=True)
    dx, dy = 2.0 / 8.0, 4.0 / 7.0

    assert phi.getValueOnGrid(*source) == pytest.approx(_cell_average(dx, dy, 0, 0), rel=2e-13)
    assert phi.getValueOnGrid(source[0] + 2, source[1] - 1) == pytest.approx(
        _cell_average(dx, dy, 2, 1), rel=2e-13
    )


def test_mode_switching_preserves_default_point_kernel():
    solver, point_phi, source = _solve()
    point_value = point_phi.getValueOnGrid(source[0] + 1, source[1])
    assert solver.getUseIntegratedGreenFunction() is False
    assert point_phi.getValueOnGrid(*source) == pytest.approx(0.0)

    solver.setUseIntegratedGreenFunction(True)
    assert solver.getUseIntegratedGreenFunction() is True
    solver.setUseIntegratedGreenFunction(False)
    assert solver.getUseIntegratedGreenFunction() is False

    rho = Grid2D(9, 8, -1.0, 1.0, -2.0, 2.0)
    phi = Grid2D(9, 8, -1.0, 1.0, -2.0, 2.0)
    rho.setValue(1.0, *source)
    solver.findPotential(rho, phi)
    assert phi.getValueOnGrid(source[0] + 1, source[1]) == pytest.approx(point_value, rel=0.0, abs=1e-15)
    with pytest.raises(TypeError):
        solver.setUseIntegratedGreenFunction()


def test_tracking_calculators_forward_integrated_green_function_mode():
    for calculator_type in (SpaceChargeCalc2p5D, SpaceChargeCalc2p5Drb, SpaceChargeCalcSliceBySlice2D):
        calculator = calculator_type(9, 8, 7)
        assert calculator.getUseIntegratedGreenFunction() is False
        calculator.setUseIntegratedGreenFunction(True)
        assert calculator.getUseIntegratedGreenFunction() is True
