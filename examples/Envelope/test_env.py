import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.lattice import AccNode
from orbit.lattice import AccLattice
from orbit.teapot import QuadTEAPOT, KickTEAPOT, TiltTEAPOT
from orbit.teapot import BendTEAPOT
from orbit.teapot import DriftTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.utils.consts import mass_proton
from orbit.envelope import Envelope
from orbit.envelope import EnvelopeTracker


def calc_bunch_cov(bunch: Bunch) -> np.ndarray:
    twiss_calc = BunchTwissAnalysis()
    twiss_calc.analyzeBunch(bunch)

    cov_matrix = np.zeros((6, 6))
    for i in range(6):
        for j in range(i + 1):
            cov_matrix[i, j] = twiss_calc.getCorrelation(j, i)
            cov_matrix[j, i] = cov_matrix[i, j]
    return cov_matrix


def make_lattice(nodes: list[AccNode]) -> AccLattice:
    lattice = TEAPOT_Lattice()
    for node in nodes:
        lattice.addNode(node)
    lattice.initialize()
    for node in lattice.getNodes():
        try:
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
        except:
            pass
    return lattice


def track_and_compare_rms(
    lattice: AccLattice,
    kin_energy: float,
    cov_matrix: np.ndarray,
    nparts: int = 100_000,
    verbose: int = 1,
) -> dict:
    """Track bunch/envelope and compare rms beam sizes.

    Args:
        lattice: Accelerator lattice.
        kin_energy: Synchronous particle kinetic energy [GeV].
        cov_matrix: 6 x 6 covariance matrix.
        nparts: Number of particles in bunch.
        verbose: Whether to print results.
    """
    cov_scale = 1e6

    data = {}
    for k1 in ["env", "bunch"]:
        data[k1] = {}
        for k2 in ["rms", "cov"]:
            data[k1][k2] = {}
            for k3 in ["env", "bunch"]:
                data[k1][k2][k3] = {}

    # Initialize bunch
    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(kin_energy)

    # Track bunch
    particles = np.random.multivariate_normal(np.zeros(6), cov_matrix, size=nparts)
    for x in particles:
        bunch.addParticle(*x)

    # Covariance matrix calculated from particles will be slightly different.
    cov_matrix = calc_bunch_cov(bunch)

    data["bunch"]["cov"]["in"] = cov_scale * calc_bunch_cov(bunch)
    lattice.trackBunch(bunch)
    data["bunch"]["cov"]["out"] = cov_scale * calc_bunch_cov(bunch)

    # Track envelope
    envelope = Envelope(sync_part=bunch.getSyncParticle(), cov_matrix=cov_matrix)
    envelope_tracker = EnvelopeTracker(lattice=lattice)

    data["env"]["cov"]["in"] = cov_scale * envelope.cov()
    envelope_tracker.track(envelope)
    data["env"]["cov"]["out"] = cov_scale * envelope.cov()

    # Compare
    for mode in ["env", "bunch"]:
        for loc in ["in", "out"]:
            data[mode]["rms"][loc] = np.sqrt(np.diag(data[mode]["cov"][loc]))

    if verbose:
        dims = ["x", "xp", "y", "yp", "z", "dE"]
        for key in ["in", "out"]:
            print(key.upper())
            for i in range(6):
                print("  rms {}:".format(dims[i]))
                print("    env:   {}".format(data["env"]["rms"][key][i]))
                print("    bunch: {}".format(data["bunch"]["rms"][key][i]))

    assert np.all(np.isclose(data["env"]["cov"]["in"], data["bunch"]["cov"]["in"]))

    atol = np.ones(6)
    atol[0:4] = 1e-3  # [mm mrad]
    atol[4] = 1e-3  # [mm]
    atol[5] = 1e-3  # [MeV]

    for i in range(6):
        x = data["env"]["rms"]["out"][i]
        y = data["bunch"]["rms"]["out"][i]
        assert np.abs(x - y) <= atol[i]


def make_default_cov_matrix(
    rms_x: float = 1e-3,
    rms_xp: float = 1e-3,
    rms_y: float = 1e-3,
    rms_yp: float = 1e-3,
    rms_z: float = 1e-3,
    rms_dE: float = 1e-5,
) -> np.ndarray:
    return np.diag(np.square([rms_x, rms_xp, rms_y, rms_yp, rms_z, rms_dE]))


def test_drift(
    kin_energy: float = 0.0025, length: float = 1.0, cov_matrix: np.ndarray = None
) -> None:
    nodes = [DriftTEAPOT(length=length)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_quad(
    kin_energy: float = 0.0025,
    length: float = 1.0,
    kq: float = 1.0,
    cov_matrix: np.ndarray = None,
) -> None:
    nodes = [QuadTEAPOT(length=length, kq=kq)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_dipole(
    kin_energy: float = 0.0025,
    length: float = 1.0,
    theta: float = 20.0,
    cov_matrix: np.ndarray = None,
) -> None:
    nodes = [BendTEAPOT(length=length, theta=np.radians(theta))]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_kick(
    kin_energy: float = 0.0025,
    length: float = 0.1,
    kx: float = 0.001,
    ky: float = 0.001,
    dE: float = 0.00001,
    cov_matrix: np.ndarray = None,
) -> None:
    nodes = [KickTEAPOT(kx=kx, ky=ky, dE=dE, length=length)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)
    
    
def test_tilt(
    kin_energy: float = 0.0025,
    angle: float = 0.25 * np.pi,
    cov_matrix: np.ndarray = None,
) -> None:
    nodes = [TiltTEAPOT(angle=angle)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


if __name__ == "__main__":
    test_drift()
    test_quad()
    test_dipole()
    test_kick()
    test_tilt()

