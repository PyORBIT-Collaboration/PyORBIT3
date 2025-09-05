from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist2D
from orbit.bunch_generators import WaterBagDist2D
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT


def make_lattice(
    length: float = 5.0, fill_fraction: float = 0.5, kq: float = 0.65
) -> AccLattice:
    drift_nodes = [
        DriftTEAPOT("drift1"),
        DriftTEAPOT("drift2"),
    ]
    quad_nodes = [
        QuadTEAPOT("qf1"),
        QuadTEAPOT("qd"),
        QuadTEAPOT("qf2"),
    ]
    for node in [drift_nodes[0], quad_nodes[1], drift_nodes[1]]:
        node.setLength(length * fill_fraction * 0.50)
    for node in [quad_nodes[0], quad_nodes[2]]:
        node.setLength(length * fill_fraction * 0.25)

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

    for node in lattice.getNodes():
        try:
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
        except:
            pass
        
    return lattice
