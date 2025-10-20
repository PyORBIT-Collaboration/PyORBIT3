from orbit.core.bunch import Bunch
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import teapot


def test_drift():
    length = 1.0
    nparts = 5
    node = teapot.DriftTEAPOT(name="name", length=length, nparts=nparts)
    assert node.getLength() == length
    assert node.getnParts() == nparts


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


def test_bunch_wrap():
    ringlength = 100.0
    node = teapot.BunchWrapTEAPOT(ringlength=ringlength)
    assert node.getParam("ring_length") == ringlength


def test_solenoid():
    B = 1.0
    length = 1.0
    node = teapot.SolenoidTEAPOT(length=length, B=B)
    assert node.getLength() == length
    assert node.getParam("B") == B


def test_multipole():
    name = "name"
    nparts = 2
    length = 1.0
    poles = [1, 2, 3]
    kls = [0.0, 1.0, 2.0]
    skews = [0, 1, 0]
    node = teapot.MultipoleTEAPOT(
        name=name,
        length=length,
        nparts=nparts,
        poles=poles,
        kls=kls,
        skews=skews,
    )
    assert node.getName() == name
    assert node.getLength() == length
    assert node.getnParts() == nparts
    assert node.getParam("kls") == kls
    assert node.getParam("skews") == skews
    assert node.getParam("poles") == poles


def test_quad():
    name = "name"
    kq = 0.5
    nparts = 10
    length = 1.0
    poles = [1, 2, 3]
    kls = [0.0, 1.0, 2.0]
    skews = [0, 1, 0]
    node = teapot.QuadTEAPOT(
        name=name,
        kq=kq,
        length=length,
        nparts=nparts,
        poles=poles,
        kls=kls,
        skews=skews,
    )
    assert node.getName() == name
    assert node.getLength() == length
    assert node.getnParts() == nparts
    assert node.getParam("kq") == kq
    assert node.getParam("kls") == kls
    assert node.getParam("skews") == skews
    assert node.getParam("poles") == poles
