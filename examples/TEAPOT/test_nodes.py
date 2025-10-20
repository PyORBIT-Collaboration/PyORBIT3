import math

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
    nparts = 5
    waveform = None

    node = teapot.SolenoidTEAPOT(length=length, nparts=nparts, B=B, waveform=waveform)
    assert node.getLength() == length
    assert node.getnParts() == nparts
    assert node.getParam("B") == B
    assert node.waveform == waveform


def test_multipole():
    name = "name"
    length = 1.0
    nparts = 2
    poles = [1, 2, 3]
    kls = [0.0, 1.0, 2.0]
    skews = [0, 1, 0]
    waveform = None

    node = teapot.MultipoleTEAPOT(
        name=name,
        length=length,
        nparts=nparts,
        poles=poles,
        kls=kls,
        skews=skews,
        waveform=waveform,
    )
    assert node.getName() == name
    assert node.getLength() == length
    assert node.getnParts() == nparts
    assert node.getParam("kls") == kls
    assert node.getParam("skews") == skews
    assert node.getParam("poles") == poles
    assert node.waveform == waveform


def test_quad():
    name = "name"
    length = 1.0
    nparts = 10
    kq = 0.5
    poles = [1, 2, 3]
    kls = [0.0, 1.0, 2.0]
    skews = [0, 1, 0]
    waveform = None

    node = teapot.QuadTEAPOT(
        name=name,
        length=length,
        nparts=nparts,
        kq=kq,
        poles=poles,
        kls=kls,
        skews=skews,
        waveform=waveform,
    )
    assert node.getName() == name
    assert node.getLength() == length
    assert node.getnParts() == nparts
    assert node.getParam("kq") == kq
    assert node.getParam("kls") == kls
    assert node.getParam("skews") == skews
    assert node.getParam("poles") == poles
    assert node.waveform == waveform


def test_bend():
    name = "name"
    length = 1.0
    nparts = 2

    poles = [1, 2, 3]
    kls = [0.0, 1.0, 2.0]
    skews = [0, 1, 0]

    ea1 = 0.0
    ea2 = 0.0
    theta = 1e-12

    node = teapot.BendTEAPOT(name=name, length=length, nparts=nparts, poles=poles, kls=kls, skews=skews, ea1=ea1, ea2=ea2, theta=theta)
    assert node.getName() == name
    assert node.getLength() == length
    assert node.getnParts() == nparts
    assert node.getParam("kls") == kls
    assert node.getParam("skews") == skews
    assert node.getParam("poles") == poles
    assert node.getParam("ea1") == ea1
    assert node.getParam("ea2") == ea2
    assert node.getParam("theta") == theta


def test_ring_rf():
    name = "name"
    harmonics = [1.0, 2.0, 3.0]
    voltages = [1e-6, 2e-6, 3e-6]
    phases = [math.pi, 0.5 * math.pi, 0.25 * math.pi]
    ringlength = 100.0

    node = teapot.RingRFTEAPOT(
        name=name,
        harmonics=harmonics,
        voltages=voltages,
        phases=phases,
        ringlength=ringlength,
    )
    assert node.getName() == name
    assert node.getParam("voltages") == voltages
    assert node.getParam("phases") == phases
    assert node.getParam("harmonics") == harmonics


def test_kick():
    name = "name"
    length = 1.0
    nparts = 4
    kx = 0.2
    ky = 0.1
    dE = 0.001
    waveform = None

    node = teapot.KickTEAPOT(
        name=name,
        length=length,
        nparts=nparts,
        kx=kx,
        ky=ky,
        dE=dE,
        waveform=waveform,
    )
    assert node.getName() == name
    assert node.getLength() == length
    assert node.getnParts() == nparts
    assert node.getParam("kx") == kx
    assert node.getParam("ky") == ky
    assert node.getParam("dE") == dE
    assert node.waveform == waveform


def test_tilt():
    name = "name"
    angle = 0.1

    node = teapot.TiltTEAPOT(angle=angle, name=name)
    assert node.getName() == name
    assert node.getTiltAngle() == angle
