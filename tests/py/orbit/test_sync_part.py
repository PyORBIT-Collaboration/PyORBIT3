from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle


def make_bunch():
    bunch = Bunch()
    bunch.mass(0.938272)
    bunch.charge(1.0)
    bunch.getSyncParticle().kinEnergy(1.0)
    return bunch


def test_get_mass_charge():
    bunch = make_bunch()
    sync_part = bunch.getSyncParticle()

    assert sync_part.mass() == bunch.mass()
    assert sync_part.charge() == bunch.charge()
