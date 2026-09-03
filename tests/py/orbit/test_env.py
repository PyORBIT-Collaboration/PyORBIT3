import numpy as np
import pytest

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.linac import MatrixRfGap
from orbit.bunch_utils import collect_bunch
from orbit.envelope import Envelope
from orbit.lattice import AccNode
from orbit.lattice import AccLattice
from orbit.py_linac.lattice import Drift
from orbit.py_linac.lattice import Quad
from orbit.py_linac.lattice import Bend
from orbit.py_linac.lattice import TiltElement
from orbit.py_linac.lattice import Solenoid
from orbit.teapot import BendTEAPOT
from orbit.teapot import ContinuousLinearFocusingTEAPOT
from orbit.teapot import DriftTEAPOT
from orbit.teapot import KickTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import SolenoidTEAPOT
from orbit.teapot import TiltTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.utils.consts import mass_proton


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
    charge: float = 1.0,
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
            for k3 in ["in", "out"]:
                data[k1][k2][k3] = {}

    # Initialize bunch
    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.charge(charge)

    sync_part = bunch.getSyncParticle()
    sync_part.kinEnergy(kin_energy)

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
    envelope = Envelope(sync_part=sync_part, cov_matrix=cov_matrix)

    data["env"]["cov"]["in"] = cov_scale * envelope.cov_matrix
    lattice.trackEnvelope(envelope)
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


def test_drift_teapot():
    node = DriftTEAPOT(length=1.0, nparts=6)
    lattice = make_lattice([node])
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix)


def test_drift_linac():
    node = Drift()
    node.setLength(1.0)
    node.setnParts(6)
    nodes = [node]

    lattice = make_lattice(nodes)
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_quad_teapot(charge: float):
    node = QuadTEAPOT(length=1.0, kq=1.0, nparts=10)
    lattice = make_lattice([node])
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix, charge=charge)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_cf_teapot(charge: float):
    node = ContinuousLinearFocusingTEAPOT(length=10.0, kq=1.0, nparts=10)
    lattice = make_lattice([node])
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix, charge=charge)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_quad_linac(charge: float):
    node = Quad()
    node.setLength(1.0)
    node.setnParts(10)
    node.setParam("dB/dr", 0.23)
    nodes = [node]

    lattice = make_lattice(nodes)
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix, charge=charge)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_bend_teapot(charge: float):
    node = BendTEAPOT(length=1.0, theta=np.radians(20.0), nparts=5)
    lattice = make_lattice([node])
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix, charge=charge)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_bend_linac(charge: float):
    node = Bend()
    node.setLength(1.0)
    node.setnParts(5)
    node.setParam("theta", np.radians(20.0))
    nodes = [node]

    lattice = make_lattice(nodes)
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix, charge=charge)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_kick_teapot(charge: float):
    node = KickTEAPOT(kx=0.001, ky=0.001, dE=0.00001, length=0.1, nparts=4)
    lattice = make_lattice([node])
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix, charge=charge)


def test_tilt_teapot():
    node = TiltTEAPOT(angle=(0.25 * np.pi))
    lattice = make_lattice([node])
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix)


def test_tilt_linac():
    node = TiltElement()
    node.setTiltAngle(0.25 * np.pi)
    lattice = make_lattice([node])
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_solenoid_teapot(charge: float):
    node = SolenoidTEAPOT(length=2.0, B=1.0, nparts=10)
    lattice = make_lattice([node])
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix, charge=charge)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_solenoid_linac(charge: float):
    node = Solenoid()
    node.setLength(2.0)
    node.setnParts(10)
    node.setParam("B", 1.0)
    nodes = [node]

    lattice = make_lattice(nodes)
    cov_matrix = make_default_cov_matrix()
    track_and_compare_rms(lattice, kin_energy=0.0025, cov_matrix=cov_matrix, charge=charge)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_rf_gap_matrix(charge: float):
    kin_energy = 0.0025
    frequency = 402.5e06
    E0TL = 0.001
    phase = 0.0

    cov_matrix = make_default_cov_matrix()

    bunch_in = Bunch()
    bunch_in.mass(mass_proton)
    bunch_in.charge(charge)
    bunch_in.getSyncParticle().kinEnergy(kin_energy)

    coords_in = np.random.multivariate_normal(np.zeros(6), cov_matrix, size=10)
    for x in coords_in:
        bunch_in.addParticle(*x)

    bunch_out_1 = Bunch()
    bunch_in.copyBunchTo(bunch_out_1)

    matrix_rf_gap = MatrixRfGap()
    matrix_rf_gap.trackBunch(bunch_out_1, frequency, E0TL, phase)

    coords_out_1 = collect_bunch(bunch_out_1)["coords"]

    from orbit.utils.matrix import get_matrix_rf_gap

    bunch_out_2 = Bunch()
    bunch_in.copyBunchTo(bunch_out_2)

    sync_part = bunch_in.getSyncParticle()
    matrix = get_matrix_rf_gap(
        sync_part,
        frequency=frequency,
        E0TL=E0TL,
        phase=phase,
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
    pass


def test_track_sublattice_no_error():
    bunch = Bunch()
    bunch.mass(mass_proton)

    sync_part = bunch.getSyncParticle()
    sync_part.kinEnergy(0.001)

    cov_matrix = np.diag(np.square([1e-3, 0, 1e-3, 0.0, 1e-3, 0.0]))
    envelope = Envelope(sync_part=sync_part, cov_matrix=cov_matrix)

    lattice = TEAPOT_Lattice()

    n = 5
    for _ in range(n):
        lattice.addNode(DriftTEAPOT(length=0.1))

    for i in range(n):
        lattice.trackEnvelope(envelope, index_start=i)
        lattice.trackEnvelope(envelope, index_stop=-i)


@pytest.mark.parametrize("charge", [1.0, -1.0])
def test_get_total_matrix(charge: float) -> None:
    node = DriftTEAPOT(length=2.0, nparts=50)
    lattice = make_lattice([node])

    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.charge(charge)

    sync_part = bunch.getSyncParticle()
    sync_part.kinEnergy(0.001)

    cov_matrix = make_default_cov_matrix()
    envelope = Envelope(sync_part=sync_part, cov_matrix=cov_matrix, intensity=1e7)

    envelope_out_a = envelope.copy()
    lattice.trackEnvelope(envelope_out_a, sc="2d")

    matrix = lattice.getEnvelopeTransferMatrix(envelope.copy(), sc="2d")
    envelope_out_b = envelope.copy()
    envelope_out_b.transform(matrix)
    assert np.all(np.isclose(envelope_out_a.cov_matrix, envelope_out_b.cov_matrix))
