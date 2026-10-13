from pathlib import Path
from pprint import pprint

import pytest

from orbit.core.bunch import Bunch
from orbit.parsers import MADX_Parser
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice


def test_sns_ring_eigtunes_coupled():
    madx_file = Path(__file__).parent / "inputs" / "sns_ring_equal_tunes.lat"

    lattice = TEAPOT_Lattice()
    lattice.readMADX(str(madx_file), "rnginjsol")
    lattice.initialize()

    for name in ["scbdsol_c13a", "scbdsol_c13b"]:
        node = lattice.getNodeForName(name)
        node.setParam("B", 0.15 / (2.0 * node.getLength()))
        
    bunch = Bunch()
    bunch.mass(0.938)
    bunch.getSyncParticle().kinEnergy(1.3)

    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
    ring_params = matrix_lattice.getRingParametersDict()

    nux = ring_params["fractional tune x"]
    nuy = ring_params["fractional tune y"]
    nu1 = ring_params["fractional tune 1"]
    nu2 = ring_params["fractional tune 2"]
    assert abs(nux - nu1) > 0.01
    assert abs(nuy - nu2) > 0.01


def test_sns_ring_eigtunes_uncoupled():
    madx_file = Path(__file__).parent / "inputs" / "sns_ring_equal_tunes.lat"

    lattice = TEAPOT_Lattice()
    lattice.readMADX(str(madx_file), "rnginjsol")
    lattice.initialize()

    bunch = Bunch()
    bunch.mass(0.938)
    bunch.getSyncParticle().kinEnergy(1.3)

    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
    ring_params = matrix_lattice.getRingParametersDict()

    nux = ring_params["fractional tune x"]
    nuy = ring_params["fractional tune y"]
    nu1 = ring_params["fractional tune 1"]
    nu2 = ring_params["fractional tune 2"]
    assert abs(nux - nu1) < 1e-8
    assert abs(nuy - nu2) < 1e-8


