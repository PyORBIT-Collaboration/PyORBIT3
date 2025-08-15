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
    d = collect_bunch(bunch)

    n_particles = bunch.getSize()

    expected_keys = {
        "x",
        "xp",
        "y",
        "yp",
        "z",
        "dE",
        "charge",
        "classical_radius",
        "mass",
        "macro_size",
        "sync_part_coords",
        "sync_part_kin_energy",
        "sync_part_momentum",
        "sync_part_beta",
        "sync_part_gamma",
        "sync_part_time",
    }

    x, xp, y, yp, z, dE = [], [], [], [], [], []
    for i in range(n_particles):
        x.append(bunch.x(i))
        xp.append(bunch.px(i))
        y.append(bunch.y(i))
        yp.append(bunch.py(i))
        z.append(bunch.z(i))
        dE.append(bunch.dE(i))

    assert set(d.keys()) == expected_keys
    assert (d["x"] == x).all()
    assert (d["xp"] == xp).all()
    assert (d["y"] == y).all()
    assert (d["yp"] == yp).all()
    assert (d["z"] == z).all()
    assert (d["dE"] == dE).all()

    assert d["charge"] == bunch.bunchAttrDouble("charge")
    assert d["classical_radius"] == bunch.bunchAttrDouble("classical_radius")
    assert d["mass"] == bunch.bunchAttrDouble("mass")
    assert d["macro_size"] == bunch.bunchAttrDouble("macro_size")

    sync_part = bunch.getSyncParticle()
    assert (d["sync_part_coords"] == sync_part.pVector()).all()
    assert d["sync_part_kin_energy"] == sync_part.kinEnergy()
    assert d["sync_part_momentum"] == sync_part.momentum()
    assert d["sync_part_beta"] == sync_part.beta()
    assert d["sync_part_gamma"] == sync_part.gamma()
    assert d["sync_part_time"] == sync_part.time()


def test_collect_empty_bunch():
    bunch = Bunch()
    d = collect_bunch(bunch)
    assert len(d) == 0


def test_collect_arbitrary_bunch_attr(bunch):
    bunch.bunchAttrDouble("arbitrary_dbl_attr", 42.0)
    bunch.bunchAttrInt("arbitrary_int_attr", 42)

    d = collect_bunch(bunch)

    assert "arbitrary_dbl_attr" in d
    assert d["arbitrary_dbl_attr"] == 42.0
    assert "arbitrary_int_attr" in d
    assert d["arbitrary_int_attr"] == 42
