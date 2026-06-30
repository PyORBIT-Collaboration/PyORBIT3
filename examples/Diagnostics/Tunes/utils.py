import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.lattice import AccLattice
from orbit.teapot import TEAPOT_Ring
from orbit.teapot import TEAPOT_MATRIX_Lattice


def get_tmat(lattice: AccLattice, bunch: Bunch) -> np.ndarray:
    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)

    M = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            M[i, j] = matrix_lattice.getOneTurnMatrix().get(i, j)
    return M


def make_lattice_sns(sol: bool = False) -> TEAPOT_Ring:
    lattice = TEAPOT_Ring()
    lattice.readMADX("inputs/sns_ring.lat", "rnginjsol")
    lattice.initialize()

    for node in lattice.getNodes():
        try:
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
        except:
            pass

        for name in ["scbdsol_c13a", "scbdsol_c13b"]:
            node = lattice.getNodeForName(name)
            B = 0.0
            if sol:
                B = 0.15 / 2.0
            node.setParam("B", B)
    return lattice