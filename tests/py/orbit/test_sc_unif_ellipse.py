import math
import random
import time

import pytest

from orbit.bunch_generators import KVDist2D
from orbit.bunch_generators import TwissContainer
from orbit.core.bunch import Bunch
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse2D
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse
from orbit.space_charge.sc2p5d import SC2p5D_AccNode
from orbit.space_charge.sc2p5d import SCUnifEllipse2D_AccNode
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.space_charge.sc2p5d import setSCUnifEllipse2DAccNodes
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import TEAPOT_Lattice


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


def _make_fodo_lattice(name):
    lattice = TEAPOT_Lattice(name)
    for index in range(2):
        nodes = (
            QuadTEAPOT(f"qf_{index}", length=0.20, kq=+1.0),
            DriftTEAPOT(f"drift_1_{index}", length=0.40),
            QuadTEAPOT(f"qd_{index}", length=0.40, kq=-1.0),
            DriftTEAPOT(f"drift_2_{index}", length=0.40),
            QuadTEAPOT(f"qf_end_{index}", length=0.20, kq=+1.0),
        )
        for node in nodes:
            node.setnParts(2)
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
            lattice.addNode(node)
    lattice.initialize()
    return lattice


def _make_kv_bunch(n_particles=16_384):
    bunch = Bunch()
    bunch.mass(0.93827231)
    bunch.charge(1.0)
    bunch.getSyncParticle().kinEnergy(1.0)

    distribution = KVDist2D(
        TwissContainer(alpha=0.0, beta=2.0, emittance=1.0e-6),
        TwissContainer(alpha=0.0, beta=2.0, emittance=1.0e-6),
    )
    random_state = random.getstate()
    random.seed(1_234_567)
    try:
        for _ in range(n_particles // 2):
            x, xp, y, yp = distribution.getCoordinates()
            z = 0.10 * (random.random() - 0.5)
            bunch.addParticle(x, xp, y, yp, z, 0.0)
            bunch.addParticle(-x, -xp, -y, -yp, -z, 0.0)
    finally:
        random.setstate(random_state)

    bunch.macroSize(2.0e11 / bunch.getSizeGlobal())
    return bunch


def _rms(bunch, coordinate):
    values = [getattr(bunch, coordinate)(index) for index in range(bunch.getSize())]
    average = sum(values) / len(values)
    return math.sqrt(sum((value - average) ** 2 for value in values) / len(values))


def _relative_particle_difference(reference, candidate, coordinate):
    reference_values = [
        getattr(reference, coordinate)(index) for index in range(reference.getSize())
    ]
    candidate_values = [
        getattr(candidate, coordinate)(index) for index in range(candidate.getSize())
    ]
    average = sum(reference_values) / len(reference_values)
    scale = math.sqrt(
        sum((value - average) ** 2 for value in reference_values)
        / len(reference_values)
    )
    difference = math.sqrt(
        sum(
            (reference_value - candidate_value) ** 2
            for reference_value, candidate_value in zip(
                reference_values, candidate_values
            )
        )
        / len(reference_values)
    )
    return difference / scale


def test_uniform_ellipse_2d_matches_grid_solver_for_kv_beam():
    initial_bunch = _make_kv_bunch()
    grid_bunch = Bunch()
    ellipse_bunch = Bunch()
    initial_bunch.copyBunchTo(grid_bunch)
    initial_bunch.copyBunchTo(ellipse_bunch)

    grid_lattice = _make_fodo_lattice("grid")
    ellipse_lattice = _make_fodo_lattice("uniform ellipse")
    path_length_min = 0.10

    grid_nodes = setSC2p5DAccNodes(
        grid_lattice,
        path_length_min,
        SpaceChargeCalc2p5D(64, 64, 1),
    )
    ellipse_nodes = setSCUnifEllipse2DAccNodes(
        ellipse_lattice,
        path_length_min,
        SpaceChargeCalcUnifEllipse2D(),
    )
    assert grid_nodes and all(isinstance(node, SC2p5D_AccNode) for node in grid_nodes)
    assert ellipse_nodes and all(
        isinstance(node, SCUnifEllipse2D_AccNode) for node in ellipse_nodes
    )
    assert len(grid_nodes) == len(ellipse_nodes)

    start = time.perf_counter()
    grid_lattice.trackBunch(grid_bunch)
    grid_time = time.perf_counter() - start
    start = time.perf_counter()
    ellipse_lattice.trackBunch(ellipse_bunch)
    ellipse_time = time.perf_counter() - start

    grid_rms = {coordinate: _rms(grid_bunch, coordinate) for coordinate in ("x", "y")}
    ellipse_rms = {
        coordinate: _rms(ellipse_bunch, coordinate) for coordinate in ("x", "y")
    }
    for coordinate in ("x", "y"):
        assert ellipse_rms[coordinate] == pytest.approx(
            grid_rms[coordinate], rel=0.02
        )

    distribution_difference = {
        coordinate: _relative_particle_difference(
            grid_bunch, ellipse_bunch, coordinate
        )
        for coordinate in ("x", "xp", "y", "yp")
    }
    assert all(value < 0.05 for value in distribution_difference.values())

    print(
        f"SC2p5D: {grid_time:.3f} s, "
        f"SCUnifEllipse2D: {ellipse_time:.3f} s, "
        f"rms: {grid_rms} vs. {ellipse_rms}, "
        f"distribution difference: {distribution_difference}"
    )
