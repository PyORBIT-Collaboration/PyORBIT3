import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.lattice import AccNode
from orbit.lattice import AccLattice
from orbit.teapot import QuadTEAPOT
from orbit.teapot import BendTEAPOT
from orbit.teapot import DriftTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.utils.consts import mass_proton

from orbit.envelope import Envelope
from orbit.envelope import EnvelopeTracker


def calc_bunch_cov(bunch: Bunch) -> np.ndarray:
    twiss_calc = BunchTwissAnalysis()
    twiss_calc.computeBunchMoments(bunch, 2, 0, 0)

    cov_matrix = np.zeros((6, 6))
    for i in range(6):
        for j in range(i + 1):
            cov_matrix[i, j] = twiss_calc.getCorrelation(j, i)
            cov_matrix[j, i] = cov_matrix[i, j]
    return cov_matrix


def track_and_compare(
    nodes: list[AccNode],
    kin_energy: float,
    cov_matrix: np.ndarray,
    nparts: int = 100_000,
) -> dict:

    cov_scale = 1e6

    data = {}
    for k1 in ["env", "bunch"]:
        data[k1] = {}
        for k2 in ["rms", "cov"]:
            data[k1][k2] = {}
            for k3 in ["env", "bunch"]:
                data[k1][k2][k3] = {}

    lattice = TEAPOT_Lattice()
    for node in nodes:
        lattice.addNode(node)
    lattice.initialize()
    for node in lattice.getNodes():
        node.setUsageFringeFieldIN(False)
        node.setUsageFringeFieldOUT(False)

    bunch = Bunch()
    bunch.mass(mass_proton)

    sync_part = bunch.getSyncParticle()
    sync_part.kinEnergy(kin_energy)

    envelope = Envelope(sync_part=sync_part, cov_matrix=cov_matrix)

    envelope_tracker = EnvelopeTracker(lattice=lattice)

    data["env"]["cov"]["in"] = cov_scale * envelope.cov()
    envelope_tracker.track(envelope)
    data["env"]["cov"]["out"] = cov_scale * envelope.cov()

    particles = np.random.multivariate_normal(np.zeros(6), cov_matrix, size=nparts)
    for x in particles:
        bunch.addParticle(*x)

    data["bunch"]["cov"]["in"] = cov_scale * calc_bunch_cov(bunch)
    lattice.trackBunch(bunch)
    data["bunch"]["cov"]["out"] = cov_scale * calc_bunch_cov(bunch)

    for mode in ["env", "bunch"]:
        for loc in ["in", "out"]:
            data[mode]["rms"][loc] = np.sqrt(np.diag(data[mode]["cov"][loc]))

    dims = ["x", "xp", "y", "yp", "z", "dE"]
    for key in ["in", "out"]:
        print(key.upper())
        for i in range(6):
            print("  rms {}:".format(dims[i]))
            print("    env:   {}".format(data["env"]["rms"][key][i]))
            print("    bunch: {}".format(data["bunch"]["rms"][key][i]))

    for key in ["in", "out"]:
        assert np.all(np.isclose(data["env"]["cov"][key], data["bunch"]["cov"][key]))


def make_default_cov_matrix(scale: float = 0.001) -> np.ndarray:
    cov_matrix = np.zeros((6, 6))
    cov_matrix[0, 0] = scale ** 2
    cov_matrix[2, 2] = scale ** 2
    return cov_matrix


def test_drift(kin_energy: float = 0.0025, length: float = 1.0, cov_matrix: np.ndarray = None):
    nodes = [
        DriftTEAPOT(length=1.0),
    ]
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare(nodes, kin_energy, cov_matrix)

def test_quad(kin_energy: float = 0.0025, length: float = 1.0, kq: float = 1.0, cov_matrix: np.ndarray = None):
    nodes = [
        QuadTEAPOT(length=length, kq=kq),
    ]
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare(nodes, kin_energy, cov_matrix)


def test_dipole(kin_energy: float = 0.0025, length: float = 1.0, theta: float = 20.0, cov_matrix: np.ndarray = None):
    nodes = [
        BendTEAPOT(length=length, theta=np.radians(theta))
    ]
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare(nodes, kin_energy, cov_matrix)


if __name__ == "__main__":
    functions = [
        # test_drift,
        # test_quad,
        test_dipole,
    ]
    for function in functions:
        function()
