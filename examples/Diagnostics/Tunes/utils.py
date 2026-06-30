import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist2D
from orbit.bunch_generators import WaterBagDist2D
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_Ring
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT


def get_tmat(lattice: AccLattice, bunch: Bunch) -> np.ndarray:
    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)

    M = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            M[i, j] = matrix_lattice.getOneTurnMatrix().get(i, j)
    return M

def make_lattice_fodo(
    length: float = 5.0,
    fill_frac: float = 0.5,
    kq: float = 0.65,
    start: str = "drift",
) -> TEAPOT_Ring:
    """Create FODO lattice.

    Args:
        length: Length of lattice [m].
        fill_frac: Fraction of lattice occupied by quadrupoles.
        kq: Quad coefficient [1/m].
        start: Whether to start in drift or quad center. {"drift", "quad"}
    """
    length_quad = length * fill_frac / 2.0
    length_drift = length * (1.0 - fill_frac) / 2.0

    if start == "quad":
        drift_nodes = [
            DriftTEAPOT("drift1"),
            DriftTEAPOT("drift2"),
        ]
        quad_nodes = [
            QuadTEAPOT("qf1"),
            QuadTEAPOT("qd"),
            QuadTEAPOT("qf2"),
        ]

        drift_nodes[0].setLength(length_drift)
        drift_nodes[1].setLength(length_drift)

        quad_nodes[0].setLength(length_quad * 0.5)
        quad_nodes[1].setLength(length_quad)
        quad_nodes[2].setLength(length_quad * 0.5)

        quad_nodes[0].setParam("kq", +kq)
        quad_nodes[1].setParam("kq", -kq)
        quad_nodes[2].setParam("kq", +kq)

        lattice = TEAPOT_Ring()
        lattice.addNode(quad_nodes[0])
        lattice.addNode(drift_nodes[0])
        lattice.addNode(quad_nodes[1])
        lattice.addNode(drift_nodes[1])
        lattice.addNode(quad_nodes[2])
        lattice.initialize()

    elif start == "drift":

        drift_nodes = [
            DriftTEAPOT("drift1"),
            DriftTEAPOT("drift2"),
        ]
        quad_nodes = [
            QuadTEAPOT("qf1"),
            QuadTEAPOT("qd"),
            QuadTEAPOT("qf2"),
        ]

        drift_nodes[0].setLength(length_drift)
        drift_nodes[1].setLength(length_drift)

        quad_nodes[0].setLength(length_quad * 0.5)
        quad_nodes[1].setLength(length_quad)
        quad_nodes[2].setLength(length_quad * 0.5)

        quad_nodes[0].setParam("kq", +kq)
        quad_nodes[1].setParam("kq", -kq)
        quad_nodes[2].setParam("kq", +kq)

        lattice = TEAPOT_Lattice()
        lattice.addNode(quad_nodes[0])
        lattice.addNode(drift_nodes[0])
        lattice.addNode(quad_nodes[1])
        lattice.addNode(drift_nodes[1])
        lattice.addNode(quad_nodes[2])
        lattice.initialize()

    else:
        raise ValueError

    for node in lattice.getNodes():
        node.setUsageFringeFieldIN(False)
        node.setUsageFringeFieldOUT(False)

    return lattice


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