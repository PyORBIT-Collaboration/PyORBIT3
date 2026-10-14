import math

import pytest

from orbit.core.spacecharge import Grid3D, PoissonSolverFFT3D, SpaceChargeCalc3D


def _antiderivative(x, y, z):
    radius = math.hypot(x, y, z)

    def log_argument(coordinate, other_a, other_b):
        if coordinate >= 0.0:
            return math.log(coordinate + radius)
        return math.log((other_a * other_a + other_b * other_b) / (radius - coordinate))

    value = 0.0
    if y and z:
        value += y * z * log_argument(x, y, z)
    if x and z:
        value += x * z * log_argument(y, x, z)
    if x and y:
        value += x * y * log_argument(z, x, y)
    if z:
        value -= z * z * math.atan(x * y / (z * radius)) / 2.0
    if y:
        value -= y * y * math.atan(x * z / (y * radius)) / 2.0
    if x:
        value -= x * x * math.atan(y * z / (x * radius)) / 2.0
    return value


def _cell_average(dx, dy, dz, x, y, z, image=0.0):
    def f(ix, iy, iz):
        return _antiderivative((ix - 0.5) * dx, (iy - 0.5) * dy, (iz - 0.5) * dz + image)

    integral = (
        f(x + 1, y + 1, z + 1) - f(x, y + 1, z + 1)
        - f(x + 1, y, z + 1) + f(x, y, z + 1)
        - f(x + 1, y + 1, z) + f(x, y + 1, z)
        + f(x + 1, y, z) - f(x, y, z)
    )
    return integral / (dx * dy * dz)


def _solve(use_integrated=False, external_bunches=0, spacing=math.inf):
    shape = (9, 8, 7)
    limits = (-1.0, 1.0, -2.0, 2.0, -3.5, 3.5)
    solver = PoissonSolverFFT3D(*shape, *limits)
    if external_bunches:
        solver.numExtBunches(external_bunches)
        solver.distBetweenBunches(spacing)
        solver.updateGeenFunction()
    solver.setUseIntegratedGreenFunction(use_integrated)

    rho = Grid3D(*shape)
    phi = Grid3D(*shape)
    for grid in (rho, phi):
        grid.setGridX(limits[0], limits[1])
        grid.setGridY(limits[2], limits[3])
        grid.setGridZ(limits[4], limits[5])
    source = (4, 4, 3)
    rho.setValue(1.0, *source)
    solver.findPotential(rho, phi)
    return solver, phi, source


def test_integrated_kernel_matches_direct_cell_integral_on_rectangular_grid():
    _, phi, source = _solve(use_integrated=True)
    dx, dy, dz = 2.0 / 8.0, 4.0 / 7.0, 7.0 / 7.0

    assert phi.getValueOnGrid(*source) == pytest.approx(_cell_average(dx, dy, dz, 0, 0, 0), rel=2e-13)
    assert phi.getValueOnGrid(source[0] + 2, source[1] - 1, source[2] + 1) == pytest.approx(
        _cell_average(dx, dy, dz, 2, 1, 1), rel=2e-13
    )


def test_mode_switching_preserves_default_point_kernel_and_python_api():
    solver, point_phi, source = _solve()
    point_value = point_phi.getValueOnGrid(source[0] + 1, source[1], source[2])
    assert solver.getUseIntegratedGreenFunction() is False
    assert point_phi.getValueOnGrid(*source) == pytest.approx(0.0)
    assert point_value == pytest.approx(4.0)

    solver.setUseIntegratedGreenFunction(True)
    assert solver.getUseIntegratedGreenFunction() is True
    solver.setUseIntegratedGreenFunction(False)
    assert solver.getUseIntegratedGreenFunction() is False
    rho = Grid3D(9, 8, 7)
    phi = Grid3D(9, 8, 7)
    for grid in (rho, phi):
        grid.setGridX(-1.0, 1.0)
        grid.setGridY(-2.0, 2.0)
        grid.setGridZ(-3.5, 3.5)
    rho.setValue(1.0, *source)
    solver.findPotential(rho, phi)
    assert phi.getValueOnGrid(source[0] + 1, source[1], source[2]) == pytest.approx(point_value, rel=0.0, abs=1e-15)

    calc = SpaceChargeCalc3D(9, 8, 7)
    assert calc.getUseIntegratedGreenFunction() is False
    calc.setUseIntegratedGreenFunction(True)
    assert calc.getUseIntegratedGreenFunction() is True
    with pytest.raises(TypeError):
        solver.setUseIntegratedGreenFunction()


def test_integrated_neighboring_bunches_sum_cell_integrals():
    spacing = 5.0
    _, phi, source = _solve(use_integrated=True, external_bunches=2, spacing=spacing)
    dx, dy, dz = 2.0 / 8.0, 4.0 / 7.0, 7.0 / 7.0
    expected = sum(_cell_average(dx, dy, dz, 0, 0, 0, image) for image in (-spacing, 0.0, spacing))
    assert math.isfinite(phi.getValueOnGrid(*source))
    assert phi.getValueOnGrid(*source) == pytest.approx(expected, rel=2e-13)
