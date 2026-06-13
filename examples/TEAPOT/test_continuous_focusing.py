import random

from orbit.core import teapot_base
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import teapot
from orbit.utils.consts import mass_proton


def test_clf_quad_teapot_base():
    """Test that CF track function gives same results as focusing quadrupole."""
    coords = []
    for i in range(100):
        x = random.gauss(mu=0.0, sigma=0.010)
        y = random.gauss(mu=0.0, sigma=0.010)
        z = random.gauss(mu=0.0, sigma=1.0)
        xp = random.gauss(mu=0.0, sigma=0.010)
        yp = random.gauss(mu=0.0, sigma=0.010)
        dE = random.gauss(mu=0.0, sigma=0.001)
        coords.append([x, xp, y, yp, z, dE])

    bunch_init = Bunch()
    bunch_init.mass(mass_proton)
    bunch_init.getSyncParticle().kinEnergy(1.000)

    bunch1 = Bunch()
    bunch2 = Bunch()
    bunch_init.copyBunchTo(bunch1)
    bunch_init.copyBunchTo(bunch2)

    for i in range(100):
        bunch1.addParticle(*coords[i])
        bunch2.addParticle(*coords[i])

    length = 1.0
    kq = 0.5

    teapot_base.continuousLinear(bunch1, length, kq)
    teapot_base.quad1(bunch2, length, kq)

    for i in range(bunch_init.getSize()):
        assert bunch1.x(i) == bunch2.x(i)
        assert bunch1.y(i) == bunch2.x(i)
        assert bunch1.z(i) == bunch2.z(i)
        assert bunch1.xp(i) == bunch2.xp(i)
        assert bunch1.yp(i) == bunch2.xp(i)
        assert bunch1.dE(i) == bunch2.dE(i)


def test_clf_quad_node():
    """Test that CF node gives same results as focusing quadrupole."""
    coords = []
    for i in range(100):
        x = random.gauss(mu=0.0, sigma=0.010)
        y = random.gauss(mu=0.0, sigma=0.010)
        z = random.gauss(mu=0.0, sigma=1.0)
        xp = random.gauss(mu=0.0, sigma=0.010)
        yp = random.gauss(mu=0.0, sigma=0.010)
        dE = random.gauss(mu=0.0, sigma=0.001)
        coords.append([x, xp, y, yp, z, dE])

    bunch_init = Bunch()
    bunch_init.mass(mass_proton)
    bunch_init.getSyncParticle().kinEnergy(1.000)

    bunch1 = Bunch()
    bunch2 = Bunch()
    bunch_init.copyBunchTo(bunch1)
    bunch_init.copyBunchTo(bunch2)

    for i in range(100):
        bunch1.addParticle(*coords[i])
        bunch2.addParticle(*coords[i])

    length = 1.0
    kq = 0.5
    nparts = 10

    node1 = teapot.ContinuousLinearFocusingTEAPOT(length=length, nparts=nparts, kq=kq)
    node2 = teapot.QuadTEAPOT(length=length, nparts=nparts, kq=kq)

    node1.trackBunch(bunch1)
    node2.trackBunch(bunch2)

    for i in range(bunch_init.getSize()):
        assert bunch1.x(i) == bunch2.x(i)
        assert bunch1.y(i) == bunch2.x(i)
        assert bunch1.z(i) == bunch2.z(i)
        assert bunch1.xp(i) == bunch2.xp(i)
        assert bunch1.yp(i) == bunch2.xp(i)
        assert bunch1.dE(i) == bunch2.dE(i)
