from orbit.core.bunch import Bunch
from orbit.bunch_generators import GaussDist3D
from orbit.bunch_utils import collect_bunch

from pytest import fixture


@fixture
def bunch():
    bunch = Bunch()
    bunch.mass(0.939294)
    bunch.charge(-1.0)
    bunch.getSyncParticle().kinEnergy(0.0025)
    gauss_dist = GaussDist3D()
    for i in range(10):
        x, xp, y, yp, z, dE = gauss_dist.getCoordinates()
        bunch.addParticle(x, xp, y, yp, z, dE)
    bunch.macroSize(10)
    return bunch


def test_collect_bunch(bunch):
    d = collect_bunch(bunch, return_memmap=False)

    n_particles = bunch.getSize()

    toplevel_keys = {"coords", "sync_part", "attributes"}

    attribute_keys = {"charge", "classical_radius", "mass", "macro_size"}
    sync_part_keys = {"coords", "kin_energy", "momentum", "beta", "gamma", "time"}

    x, xp, y, yp, z, dE = [], [], [], [], [], []
    for i in range(n_particles):
        x.append(bunch.x(i))
        xp.append(bunch.px(i))
        y.append(bunch.y(i))
        yp.append(bunch.py(i))
        z.append(bunch.z(i))
        dE.append(bunch.dE(i))

    assert set(d.keys()) == toplevel_keys
    assert set(d["sync_part"].keys()) == sync_part_keys
    assert set(d["attributes"].keys()) == attribute_keys

    assert d["coords"].shape == (n_particles, 6)

    assert d["attributes"]["charge"] == bunch.bunchAttrDouble("charge")
    assert d["attributes"]["classical_radius"] == bunch.bunchAttrDouble(
        "classical_radius"
    )
    assert d["attributes"]["mass"] == bunch.bunchAttrDouble("mass")
    assert d["attributes"]["macro_size"] == bunch.bunchAttrDouble("macro_size")

    sync_part = bunch.getSyncParticle()
    assert (d["sync_part"]["coords"] == sync_part.pVector()).all()
    assert d["sync_part"]["kin_energy"] == sync_part.kinEnergy()
    assert d["sync_part"]["momentum"] == sync_part.momentum()
    assert d["sync_part"]["beta"] == sync_part.beta()
    assert d["sync_part"]["gamma"] == sync_part.gamma()
    assert d["sync_part"]["time"] == sync_part.time()


def test_collect_empty_bunch():
    bunch = Bunch()
    d = collect_bunch(bunch)
    assert d is None


def test_collect_arbitrary_bunch_attr(bunch):
    bunch.bunchAttrDouble("arbitrary_dbl_attr", 42.0)
    bunch.bunchAttrInt("arbitrary_int_attr", 42)

    d = collect_bunch(bunch)

    assert "arbitrary_dbl_attr" in d["attributes"]
    assert d["attributes"]["arbitrary_dbl_attr"] == 42.0
    assert "arbitrary_int_attr" in d["attributes"]
    assert d["attributes"]["arbitrary_int_attr"] == 42
