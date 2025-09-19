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
    length: float = 5.0,
    fill_frac: float = 0.5,
    kq: float = 0.65,
    start: str = "drift",
    sol: float = 0.0,
) -> AccLattice:
    """Create FODO lattice.
    
    Args:
        length: Length of lattice [m].
        fill_frac: Fraction of lattice occupied by quadrupoles.
         
    """
    if start == "drift":
        return make_lattice_drift_start(length=length, fill_frac=fill_frac, kq=kq)
    elif start == "quad":
        return make_lattice_quad_start(length=length, fill_frac=fill_frac, kq=kq)
    else:
        raise ValueError
    
    
def make_lattice_quad_start(
    length: float = 5.0, fill_frac: float = 0.5, kq: float = 0.65
) -> AccLattice:
    
    length_quad = length * fill_frac / 2.0
    length_drift = length * (1.0 - fill_frac) / 2.0

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

    for node in lattice.getNodes():
        try:
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
        except:
            pass
        
    return lattice


def make_lattice_drift_start(
    length: float = 5.0, fill_frac: float = 0.5, kq: float = 0.65
) -> AccLattice:

    length_quad = length * fill_frac / 2.0
    length_drift = length * (1.0 - fill_frac) / 2.0

    drift_nodes = [
        DriftTEAPOT("drift1"),
        DriftTEAPOT("drift2"),
        DriftTEAPOT("drift3"),
    ]
    quad_nodes = [
        QuadTEAPOT("qf"),
        QuadTEAPOT("qd"),
    ]

    
    drift_nodes[0].setLength(length_drift * 0.5)
    drift_nodes[1].setLength(length_drift)
    drift_nodes[2].setLength(length_drift * 0.5)

    quad_nodes[0].setLength(length_quad)
    quad_nodes[1].setLength(length_quad)
    
    quad_nodes[0].setParam("kq", +kq)
    quad_nodes[1].setParam("kq", -kq)

    lattice = TEAPOT_Lattice()
    lattice.addNode(drift_nodes[0])
    lattice.addNode(quad_nodes[0])
    lattice.addNode(drift_nodes[1])
    lattice.addNode(quad_nodes[1])
    lattice.addNode(drift_nodes[2])
    lattice.initialize()

    for node in lattice.getNodes():
        try:
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
        except:
            pass
        
    return lattice
