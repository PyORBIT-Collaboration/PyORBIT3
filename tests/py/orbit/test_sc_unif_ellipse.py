import math

import pytest

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse2D
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse


def test_uniform_ellipsoid_calculator_names():
    assert SpaceChargeCalcUnifEllipse2D().getNEllipses() == 1
    assert SpaceChargeCalcUnifEllipse().getNEllipses() == 1


def test_uniform_ellipse_2d_kick_matches_envelope_formula():
    bunch = Bunch()
    bunch.mass(0.93827231)
    bunch.charge(1.0)
    bunch.getSyncParticle().kinEnergy(1.0)

    radius_x = 0.004
    radius_y = 0.002
    bunch_length = 0.12
    z_rms = bunch_length / math.sqrt(12.0)
    phi = 0.37
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    local_coordinates = (
        (+radius_x / math.sqrt(2.0), 0.0, +z_rms),
        (-radius_x / math.sqrt(2.0), 0.0, -z_rms),
        (0.0, +radius_y / math.sqrt(2.0), +z_rms),
        (0.0, -radius_y / math.sqrt(2.0), -z_rms),
    )
    coordinates = [
        (cos_phi * x + sin_phi * y, -sin_phi * x + cos_phi * y, z)
        for x, y, z in local_coordinates
    ]
    for x, y, z in coordinates:
        bunch.addParticle(x, 0.0, y, 0.0, z, 0.0)

    macro_size = 2.5e8
    bunch.macroSize(macro_size)
    intensity = macro_size * len(coordinates)
    length = 0.25

    beta = bunch.getSyncParticle().beta()
    gamma = bunch.getSyncParticle().gamma()
    perveance = (
        2.0
        * intensity
        * bunch.classicalRadius()
        / (beta**2 * gamma**3 * bunch_length)
    )
    kappa_factor = 2.0 * perveance / (radius_x + radius_y)

    calculator = SpaceChargeCalcUnifEllipse2D()
    calculator.trackBunch(bunch, length)

    for index, (x, y, _) in enumerate(local_coordinates):
        xp = kappa_factor * length * x / radius_x
        yp = kappa_factor * length * y / radius_y
        expected_xp = cos_phi * xp + sin_phi * yp
        expected_yp = -sin_phi * xp + cos_phi * yp
        assert bunch.xp(index) == pytest.approx(expected_xp)
        assert bunch.yp(index) == pytest.approx(expected_yp)
    assert all(bunch.dE(index) == 0.0 for index in range(bunch.getSize()))

    x = 2.0 * radius_x
    x_lab = cos_phi * x
    y_lab = -sin_phi * x
    ex, ey, ez = calculator.calculateField(x_lab, y_lab, 0.0)
    lambda_value = x**2 - radius_x**2
    radius_y_lambda = math.sqrt(radius_y**2 + lambda_value)
    ex_local = 2.0 * (intensity / bunch_length) / (x + radius_y_lambda)
    assert ex == pytest.approx(cos_phi * ex_local)
    assert ey == pytest.approx(-sin_phi * ex_local)
    assert ez == 0.0
