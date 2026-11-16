from pathlib import Path

from orbit.core.bunch import Bunch
from orbit.parsers import MADX_Parser
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton


def matrix_is_uncoupled(matrix_lattice: TEAPOT_MATRIX_Lattice, atol: float = 1e-12):
    for i in range(4):
        for j in range(4):
            matrix_element = matrix_lattice.oneTurnMatrix.get(i, j)
            if (i >= 2 and j < 2) or (j >= 2 and i < 2):
                if abs(matrix_element) > atol:
                    return False
    return True


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
    assert not matrix_is_uncoupled(matrix_lattice)

    matrix_params = matrix_lattice.getRingParametersDict()
    nux = matrix_params["fractional tune x"]
    nuy = matrix_params["fractional tune y"]
    nu1 = matrix_params["fractional tune 1"]
    nu2 = matrix_params["fractional tune 2"]
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
    assert matrix_is_uncoupled(matrix_lattice)

    matrix_params = matrix_lattice.getRingParametersDict()
    nux = matrix_params["fractional tune x"]
    nuy = matrix_params["fractional tune y"]
    nu1 = matrix_params["fractional tune 1"]
    nu2 = matrix_params["fractional tune 2"]
    assert abs(nux - nu1) < 1e-8
    assert abs(nuy - nu2) < 1e-8


def test_frac_tunes_match_matrix_params():
    madx_file = Path(__file__).parent / "inputs" / "sns_ring.lat"

    lattice = TEAPOT_Lattice()
    lattice.readMADX(str(madx_file), "rnginjsol")
    lattice.initialize()

    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(1.3)

    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
    matrix_params = matrix_lattice.getRingParametersDict()
    tune_x, tune_y = matrix_lattice.getTunes()
    assert matrix_is_uncoupled(matrix_lattice)
    assert abs(tune_x % 1.0 - matrix_params["fractional tune x"]) < 1e-8
    assert abs(tune_y % 1.0 - matrix_params["fractional tune y"]) < 1e-8
