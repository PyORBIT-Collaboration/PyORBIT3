import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.linac import MatrixRfGap
from orbit.bunch_utils import collect_bunch
from orbit.lattice import AccNode
from orbit.lattice import AccLattice
from orbit.py_linac.lattice import Drift
from orbit.py_linac.lattice import Quad
from orbit.py_linac.lattice import Bend
from orbit.py_linac.lattice import TiltElement
from orbit.py_linac.lattice import Solenoid
from orbit.teapot import BendTEAPOT
from orbit.teapot import DriftTEAPOT
from orbit.teapot import KickTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import SolenoidTEAPOT
from orbit.teapot import TiltTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.utils.consts import mass_proton
from orbit.envelope import Envelope
from orbit.envelope import EnvelopeTracker


def get_lorentz_factors(kin_energy: float, mass: float) -> tuple[float, float]:
    gamma = 1.0 + kin_energy / mass
    beta = np.sqrt(1.0 - (1.0 / gamma) ** 2)
    return (gamma, beta)


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
    envelope = Envelope(bunch=bunch, cov_matrix=cov_matrix)
    envelope_tracker = EnvelopeTracker(lattice=lattice)

    data["env"]["cov"]["in"] = cov_scale * envelope.cov_matrix
    envelope_tracker.track(envelope)
    data["env"]["cov"]["out"] = cov_scale * envelope.cov_matrix

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


def test_drift_teapot(
    kin_energy: float = 0.0025,
    length: float = 1.0,
    cov_matrix: np.ndarray = None,
    nparts: int = 6,
) -> None:
    nodes = [DriftTEAPOT(length=length, nparts=nparts)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_drift_linac(
    kin_energy: float = 0.0025,
    length: float = 1.0,
    cov_matrix: np.ndarray = None,
    nparts: int = 6,
) -> None:
    node = Drift()
    node.setLength(length)
    node.setnParts(nparts)
    nodes = [node]

    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_quad_teapot(
    kin_energy: float = 0.0025,
    length: float = 1.0,
    kq: float = 1.0,
    cov_matrix: np.ndarray = None,
    nparts: int = 10,
) -> None:
    nodes = [QuadTEAPOT(length=length, kq=kq, nparts=nparts)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_quad_linac(
    kin_energy: float = 0.0025,
    length: float = 1.0,
    field_grad: float = 0.23,
    cov_matrix: np.ndarray = None,
    nparts: int = 10,
) -> None:
    node = Quad()
    node.setLength(length)
    node.setnParts(nparts)
    node.setParam("dB/dr", field_grad)
    nodes = [node]

    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_bend_teapot(
    kin_energy: float = 0.0025,
    length: float = 1.0,
    theta: float = 20.0,
    cov_matrix: np.ndarray = None,
    nparts: int = 2,
) -> None:
    nodes = [BendTEAPOT(length=length, theta=np.radians(theta), nparts=nparts)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_bend_linac(
    kin_energy: float = 0.0025,
    length: float = 1.0,
    theta: float = 20.0,
    cov_matrix: np.ndarray = None,
    nparts: int = 2,
) -> None:
    node = Bend()
    node.setLength(length)
    node.setnParts(nparts)
    node.setParam("theta", np.radians(theta))
    nodes = [node]

    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_kick_teapot(
    kin_energy: float = 0.0025,
    length: float = 0.1,
    kx: float = 0.001,
    ky: float = 0.001,
    dE: float = 0.00001,
    cov_matrix: np.ndarray = None,
    nparts: int = 4,
) -> None:
    nodes = [KickTEAPOT(kx=kx, ky=ky, dE=dE, length=length, nparts=nparts)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_tilt_teapot(
    kin_energy: float = 0.0025,
    angle: float = 0.25 * np.pi,
    cov_matrix: np.ndarray = None,
) -> None:
    nodes = [TiltTEAPOT(angle=angle)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_tilt_linac(
    kin_energy: float = 0.0025,
    angle: float = 0.25 * np.pi,
    cov_matrix: np.ndarray = None,
) -> None:
    node = TiltElement()
    node.setTiltAngle(angle)
    nodes = [node]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_solenoid_teapot(
    kin_energy: float = 0.0025,
    length: float = 2.0,
    B: float = 1.0,
    cov_matrix: np.ndarray = None,
    nparts: int = 10,
) -> None:
    nodes = [SolenoidTEAPOT(length=length, B=B, nparts=nparts)]
    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_solenoid_linac(
    kin_energy: float = 0.0025,
    length: float = 2.0,
    B: float = 1.0,
    cov_matrix: np.ndarray = None,
    nparts: int = 10,
) -> None:
    node = Solenoid()
    node.setLength(length)
    node.setnParts(nparts)
    node.setParam("B", B)
    nodes = [node]

    lattice = make_lattice(nodes)
    if cov_matrix is None:
        cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy, cov_matrix)


def test_rf_gap_matrix(
    kin_energy: float = 0.0025,
    frequency: float = 402.5e06,
    E0TL: float = 0.001,
    phase: float = 0.0,
):
    # Just tests matrix against MatrixRFGap. Node not implemented yet.

    cov_matrix = make_default_cov_matrix()

    bunch_in = Bunch()
    bunch_in.mass(mass_proton)
    bunch_in.getSyncParticle().kinEnergy(kin_energy)

    coords_in = np.random.multivariate_normal(np.zeros(6), cov_matrix, size=10)
    for x in coords_in:
        bunch_in.addParticle(*x)

    bunch_out_1 = Bunch()
    bunch_in.copyBunchTo(bunch_out_1)

    matrix_rf_gap = MatrixRfGap()
    matrix_rf_gap.trackBunch(bunch_out_1, frequency, E0TL, phase)

    coords_out_1 = collect_bunch(bunch_out_1)["coords"]

    from orbit.matrix_lattice.analytic import rf_gap_matrix

    bunch_out_2 = Bunch()
    bunch_in.copyBunchTo(bunch_out_2)

    matrix = rf_gap_matrix(
        frequency=frequency,
        E0TL=E0TL,
        phase=phase,
        sync_part=bunch_out_2.getSyncParticle(),
    )
    coords_in = np.column_stack([coords_in, np.ones(coords_in.shape[0])])
    coords_out_2 = np.matmul(coords_in, matrix.T)
    coords_out_2 = coords_out_2[:, :-1]
    assert np.allclose(coords_out_1, coords_out_2)


def test_sc_3d_cold_expansion():
    # This should test expansion of cold uniform-density sphere
    # (in rest frame). We can calculate the time to expand to
    # twice the initial size. (See examples from A. Shishlo or
    # from the ImpactX repo.)
    raise NotImplementedError


if __name__ == "__main__":
    for kin_energy in [0.0025, 1.0, 10.0]:
        test_drift_teapot(kin_energy=kin_energy)
        test_quad_teapot(kin_energy=kin_energy)
        test_bend_teapot(kin_energy=kin_energy)
        test_tilt_teapot(kin_energy=kin_energy)
        test_solenoid_teapot(kin_energy=kin_energy)
        test_kick_teapot(kin_energy=kin_energy)

        test_drift_linac(kin_energy=kin_energy)
        test_quad_linac(kin_energy=kin_energy)
        test_bend_linac(kin_energy=kin_energy)
        test_tilt_linac(kin_energy=kin_energy)
        test_solenoid_linac(kin_energy=kin_energy)

        test_rf_gap_matrix(kin_energy=kin_energy)

