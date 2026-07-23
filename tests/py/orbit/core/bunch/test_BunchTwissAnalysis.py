import pytest

from orbit.core.bunch import Bunch, BunchTwissAnalysis

import numpy as np

@pytest.fixture
def gaus_dist():
    return np.random.normal(loc=0, scale=1e-3, size=(100, 6))

@pytest.fixture
def gaus_bunch(gaus_dist):
    b = Bunch()
    b.getSyncParticle().kinEnergy(2.5e-3)
    b.mass(0.939294)
    b.charge(-1.0)
    for coords in gaus_dist:
        b.addParticle(*coords)
    return b

@pytest.fixture
def twiss_analysis():
    return BunchTwissAnalysis()

def test_computeBunchMoments(gaus_bunch, gaus_dist, twiss_analysis):
    twiss_analysis.computeBunchMoments(gaus_bunch)

    mean = gaus_dist.mean(axis=0)
    cov = np.cov(gaus_dist.T, bias=True)

    for i in range(6):
        assert twiss_analysis.getAverage(i) == pytest.approx(mean[i])
        for j in range(6):
            assert twiss_analysis.getCorrelation(i, j) == pytest.approx(cov[i, j])

def test_computeBunchMoments_FirstOrder(gaus_bunch, gaus_dist, twiss_analysis):
    twiss_analysis.computeBunchMoments(gaus_bunch, 3)

    mean = gaus_dist.mean(axis=0)
    cov = np.cov(gaus_dist.T, bias=True)

    assert twiss_analysis.getBunchMoment(1, 0) == pytest.approx(mean[0])
    assert twiss_analysis.getBunchMoment(0, 1) == pytest.approx(mean[2])
    assert twiss_analysis.getBunchMoment(1, 1) == pytest.approx(cov[0, 2])

def test_computeBunchMoments_SecondOrder(gaus_bunch, gaus_dist, twiss_analysis):
    twiss_analysis.computeBunchMoments(gaus_bunch, 3)

    var_x = np.var(gaus_dist[:, 0])
    var_y = np.var(gaus_dist[:, 2])

    assert twiss_analysis.getBunchMoment(2, 0) == pytest.approx(var_x)
    assert twiss_analysis.getBunchMoment(0, 2) == pytest.approx(var_y)

