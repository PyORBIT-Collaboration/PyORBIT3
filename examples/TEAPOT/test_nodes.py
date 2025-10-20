from orbit.core.bunch import Bunch
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import teapot


def test_drift():
    length = 1.0
    node = teapot.DriftTEAPOT(name="name", length=length)
    assert node.getLength() == length


def test_aperture():
    shape = 1
    dim = (1, 2)
    node = teapot.ApertureTEAPOT(name="name", shape=shape, dim=dim)
    assert node.getParam("apertype") == shape
    assert node.getParam("aperture") == dim


def test_monitor():
    name = "name"
    node = teapot.MonitorTEAPOT(name=name)
    assert node.getName() == name
