from orbit.core.bunch import Bunch
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import teapot


def test_drift():
    length = 1.0
    node = teapot.DriftTEAPOT(name="name", length=length)
    assert node.getLength() == length
