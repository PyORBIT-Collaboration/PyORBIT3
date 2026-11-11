import numpy as np
import pytest

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.linac import MatrixRfGap
from orbit.bunch_utils import collect_bunch
from orbit.envelope import Envelope
from orbit.envelope import EnvelopeAnalysis
from orbit.envelope import fit_transfer_matrix
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.lattice import AccNodeBunchTracker
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
from orbit.teapot import MultipoleTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import SolenoidTEAPOT
from orbit.teapot import TiltTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.utils.consts import mass_proton


class EnvelopeMonitor:
    def __init__(self):
        self.history = {}

    def __call__(self, params):
        envelope, s = params["envelope"], params["path_length"]
        for key, value in {"s": s, **EnvelopeAnalysis(envelope).to_dict()}.items():
            self.history.setdefault(key, []).append(value)


def monitor_actions(monitor):
    if monitor is None:
        return None
    actions = AccActionsContainer()
    actions.addAction(monitor, AccActionsContainer.ENTRANCE)
    actions.addAction(monitor, AccActionsContainer.EXIT)
    return actions


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
def test_bend_teapot_edge_matrix(charge: float):
    node = BendTEAPOT(
        length=1.0,
        theta=np.radians(20.0),
        nparts=5,
        ea1=np.radians(8.0),
        ea2=np.radians(-3.0),
    )
    lattice = make_lattice([node])

    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.charge(charge)
    bunch.getSyncParticle().kinEnergy(1.3)

    matrix_fit = fit_transfer_matrix(lattice, bunch)
    envelope = Envelope(sync_part=bunch.getSyncParticle())
    matrix_analytic = lattice.getEnvelopeTransferMatrix(envelope)
    assert np.allclose(matrix_analytic, matrix_fit, atol=1.0e-10)


def test_fitted_envelope_matrices():
    node = MultipoleTEAPOT(
        length=0.4,
        nparts=3,
        poles=[2],
        kls=[0.2],
        skews=[0],
    )
    lattice = make_lattice([node])

    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(1.3)

    matrix_expected = fit_transfer_matrix(lattice, bunch)
    envelope = Envelope(sync_part=bunch.getSyncParticle())
    matrix_actual = lattice.getEnvelopeTransferMatrix(envelope, fit=True)
    assert np.allclose(matrix_actual, matrix_expected, atol=1.0e-12)


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


def test_track_static_sublattice():
    lattice = make_lattice(
        [
            DriftTEAPOT(length=0.1),
            DriftTEAPOT(length=0.2),
            DriftTEAPOT(length=0.3),
        ]
    )

    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(0.001)
    envelope = Envelope(
        sync_part=bunch.getSyncParticle(),
        cov_matrix=make_default_cov_matrix(),
    )

    for index_start, index_stop in [(0, 0), (1, 2)]:
        expected = envelope.copy()
        lattice.trackEnvelope(
            expected,
            index_start=index_start,
            index_stop=index_stop,
            fit=False,
        )

        actual = envelope.copy()
        lattice.trackEnvelope(
            actual,
            index_start=index_start,
            index_stop=index_stop,
            fit=False,
            static=True,
        )

        assert np.allclose(actual.cov_matrix, expected.cov_matrix)
        assert np.allclose(actual.centroid, expected.centroid)


def test_track_static_history():
    lattice = make_lattice(
        [
            DriftTEAPOT(length=0.1, nparts=1),
            DriftTEAPOT(length=0.2, nparts=2),
        ]
    )

    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(0.001)
    envelope = Envelope(
        sync_part=bunch.getSyncParticle(),
        cov_matrix=make_default_cov_matrix(),
    )

    expected_monitor = EnvelopeMonitor()
    actual_monitor = EnvelopeMonitor()
    lattice.trackEnvelope(envelope.copy(), actionContainer=monitor_actions(expected_monitor), fit=False)
    lattice.trackEnvelope(
        envelope.copy(),
        actionContainer=monitor_actions(actual_monitor),
        fit=False,
        static=True,
    )

    assert actual_monitor.history.keys() == expected_monitor.history.keys()
    for key in actual_monitor.history:
        assert np.allclose(actual_monitor.history[key], expected_monitor.history[key])


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


def test_envelope_copy_reference_state():
    import gc

    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.charge(-1.0)
    sync = bunch.getSyncParticle()
    sync.rVector((1.0, 2.0, 3.0))
    sync.pVector((0.1, 0.2, 1.3))
    sync.nxVector((1.0, 0.0, 0.0))
    sync.time(2.5)
    original = Envelope(bunch=bunch)
    original.centroid = np.arange(6, dtype=float)
    original.cov_matrix = make_default_cov_matrix()
    original.intensity = 1234.0
    original.rms_bunch_length_factor = 4.0
    copied = original.copy()

    for attr in ("mass", "charge", "kinEnergy", "time", "rVector", "pVector", "nxVector", "nyVector"):
        assert np.allclose(getattr(copied.sync_part, attr)(), getattr(sync, attr)())
    assert copied.rms_bunch_length_factor == 4.0
    assert np.array_equal(copied.centroid, original.centroid)
    assert np.array_equal(copied.cov_matrix, original.cov_matrix)
    assert copied.intensity == original.intensity
    assert copied.classical_radius == original.classical_radius
    copied.centroid[0] = 9.0
    copied.cov_matrix[0, 0] = 9.0
    copied.sync_part.time(7.0)
    copied.sync_part.kinEnergy(0.5)
    assert original.centroid[0] == 0.0
    assert original.cov_matrix[0, 0] != 9.0
    assert sync.time() == 2.5
    assert sync.kinEnergy() != 0.5
    del original, bunch, sync
    gc.collect()
    assert copied.sync_part.kinEnergy() == pytest.approx(0.5)
    assert copied.copy().sync_part.time() == 7.0


@pytest.mark.parametrize("bend_class", [Bend, BendTEAPOT])
@pytest.mark.parametrize("edges", [(0.1, 0.0), (0.0, -0.1), (0.1, -0.1)])
@pytest.mark.parametrize("fringe", [False, True])
def test_edge_bends_always_fit(bend_class, edges, fringe, monkeypatch):
    node = bend_class()
    node.setLength(1.0)
    node.setnParts(5)
    node.setParam("theta", 0.2)
    node.setParam("ea1", edges[0])
    node.setParam("ea2", edges[1])
    lattice = make_lattice([DriftTEAPOT(length=0.2), node, DriftTEAPOT(length=0.3)])
    node.setUsageFringeFieldIN(fringe)
    node.setUsageFringeFieldOUT(fringe)
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    envelope = Envelope(sync_part=bunch.getSyncParticle())
    expected = fit_transfer_matrix(lattice, bunch)
    expected_per_node = lattice.getEnvelopeTransferMatrix(envelope, fit=True)

    def unexpected_analytic(*args, **kwargs):
        pytest.fail("Edge bend should use fitting even with fit=False")

    monkeypatch.setattr(node, "getMatrix", unexpected_analytic)
    actual = lattice.getEnvelopeTransferMatrix(envelope, fit=False)
    assert np.allclose(actual, expected_per_node, atol=1e-12)
    # Finite probe amplitudes give slightly different linearizations when
    # fitting each element separately versus fitting the full nonlinear map.
    assert np.allclose(actual, expected, atol=3e-8)


class UnsupportedEnvelopeDrift(DriftTEAPOT):
    def getMatrix(self, sync_part, part_index=-1):
        raise NotImplementedError("No analytic matrix")


@pytest.mark.parametrize("static", [False, True])
@pytest.mark.parametrize("history", [False, True])
def test_envelope_fit_fallback(static, history):
    lattice = make_lattice([
        DriftTEAPOT(length=0.1),
        UnsupportedEnvelopeDrift(length=0.2, nparts=2),
        DriftTEAPOT(length=0.3),
    ])
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    envelope = Envelope(sync_part=bunch.getSyncParticle())
    with pytest.raises(NotImplementedError, match="No analytic matrix"):
        lattice.trackEnvelope(
            envelope.copy(), fit=False, fallback=False, static=static,
            actionContainer=monitor_actions(EnvelopeMonitor() if history else None),
        )
    expected = envelope.copy()
    expected_monitor = EnvelopeMonitor() if history else None
    assert lattice.trackEnvelope(expected, fit=True, actionContainer=monitor_actions(expected_monitor)) is None
    for _ in range(2):  # Also exercise reuse of the static cache.
        actual = envelope.copy()
        actual_monitor = EnvelopeMonitor() if history else None
        result = lattice.trackEnvelope(
            actual, fit=False, fallback=True, static=static, actionContainer=monitor_actions(actual_monitor),
        )
        assert np.allclose(actual.cov_matrix, expected.cov_matrix)
        assert actual.sync_part.time() == pytest.approx(expected.sync_part.time())
        assert result is None
        if history:
            for key in expected_monitor.history:
                assert np.allclose(actual_monitor.history[key], expected_monitor.history[key])
    assert np.allclose(
        lattice.getEnvelopeTransferMatrix(envelope, fit=False, fallback=True),
        lattice.getEnvelopeTransferMatrix(envelope, fit=True),
    )


class AcceleratingEnvelopeNode(AccNode):
    def getMatrix(self, sync_part, part_index=-1):
        sync_part.kinEnergy(sync_part.kinEnergy() + 0.1)
        sync_part.time(sync_part.time() + 1.0)
        return np.identity(7)

    def track(self, params):
        sync = params["bunch"].getSyncParticle()
        sync.kinEnergy(sync.kinEnergy() + 0.1)
        sync.time(sync.time() + 1.0)


@pytest.mark.parametrize("sc", [None, "2d", "3d"])
@pytest.mark.parametrize("history", [False, True])
@pytest.mark.parametrize("fit", [False, True])
def test_static_envelope_sync_updates(sc, history, fit):
    lattice = make_lattice([
        AcceleratingEnvelopeNode(), UnsupportedEnvelopeDrift(length=0.1),
    ])
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    envelope = Envelope(sync_part=bunch.getSyncParticle(), intensity=1e7)
    expected = envelope.copy()
    expected_monitor = EnvelopeMonitor() if history else None
    assert lattice.trackEnvelope(expected, sc=sc, actionContainer=monitor_actions(expected_monitor)) is None
    for _ in range(2):
        actual = envelope.copy()
        actual_monitor = EnvelopeMonitor() if history else None
        result = lattice.trackEnvelope(
            actual, static=True, sc=sc, actionContainer=monitor_actions(actual_monitor), fit=fit, fallback=True,
        )
        assert np.allclose(actual.cov_matrix, expected.cov_matrix)
        assert actual.kin_energy == pytest.approx(expected.kin_energy)
        assert actual.sync_part.time() == pytest.approx(expected.sync_part.time())
        assert result is None
        if history:
            for key in expected_monitor.history:
                assert np.allclose(actual_monitor.history[key], expected_monitor.history[key])
    assert envelope.kin_energy == pytest.approx(1.3)
    lattice.getEnvelopeTransferMatrix(envelope, sc=sc)
    assert envelope.kin_energy == pytest.approx(1.3)
    assert envelope.sync_part.time() == 0.0


def test_envelope_fallback_only_handles_missing_matrices(monkeypatch):
    node = DriftTEAPOT(length=0.1)
    lattice = make_lattice([node])
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    envelope = Envelope(sync_part=bunch.getSyncParticle())

    def no_matrix(*args, **kwargs):
        return None

    monkeypatch.setattr(node, "getMatrix", no_matrix)
    unchanged = envelope.copy()
    lattice.trackEnvelope(unchanged, fit=False, fallback=True)
    assert np.array_equal(unchanged.cov_matrix, envelope.cov_matrix)
    assert unchanged.sync_part.time() == envelope.sync_part.time()

    def broken_matrix(*args, **kwargs):
        raise RuntimeError("Invalid node parameters")

    monkeypatch.setattr(node, "getMatrix", broken_matrix)
    with pytest.raises(RuntimeError, match="Invalid node parameters"):
        lattice.trackEnvelope(envelope, fit=False, fallback=True)


def test_envelope_analysis():
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    cov = np.diag([1., 4., 9., 16., 25., 36.])
    envelope = Envelope(sync_part=bunch.getSyncParticle(), cov_matrix=cov)
    analysis = EnvelopeAnalysis(envelope)
    assert np.allclose(analysis.rms(), [1., 3., 5.])
    assert analysis.rms(1) == 3.0
    assert np.allclose(analysis.emittances(), [2., 12., 30.])
    assert np.allclose(analysis.eigenemittances(), [2., 12., 30.])
    expected_normalized = np.array([2., 12., 30.]) * [
        envelope.beta * envelope.gamma, envelope.beta * envelope.gamma, 1 / envelope.beta,
    ]
    assert np.allclose(analysis.emittances(normalized=True), expected_normalized)
    data = analysis.to_dict()
    assert data["rms_z"] == 5.0
    assert data["eps_y"] == 12.0
    assert data["eps_z_n"] == pytest.approx(expected_normalized[2])
    assert data["eps_3"] == 30.0
    data["cov"][0, 0] = 99.0
    data["mean"][0] = 99.0
    assert envelope.cov_matrix[0, 0] == 1.0
    assert envelope.centroid[0] == 0.0
    envelope.cov_matrix = np.eye(6) * 4.0
    assert np.allclose(analysis.rms(), [2., 2., 2.])


def test_envelope_analysis_coupled_and_cold_modes():
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    rotation = np.eye(6)
    angle = 0.4
    for i, j in [(0, 2), (1, 3)]:
        rotation[i, i] = rotation[j, j] = np.cos(angle)
        rotation[i, j] = np.sin(angle)
        rotation[j, i] = -np.sin(angle)
    envelope = Envelope(
        sync_part=bunch.getSyncParticle(),
        cov_matrix=rotation @ np.diag([1., 4., 9., 16., 0., 0.]) @ rotation.T,
    )
    analysis = EnvelopeAnalysis(envelope)
    assert np.allclose(analysis.eigenemittances(), [0., 2., 12.])
    envelope.cov_matrix = np.zeros((6, 6))
    assert np.array_equal(analysis.eigenemittances(), np.zeros(3))
    assert analysis.to_dict()["eps_3"] == 0.0


class ActionTestNode(AccNodeBunchTracker):
    """Affine transport with observable centroid and reference-time changes."""

    def track(self, params):
        bunch = params["bunch"]
        length = self.getLength(self.getActivePartIndex())
        for i in range(bunch.getSize()):
            bunch.x(i, bunch.x(i) + length)
        sync = bunch.getSyncParticle()
        sync.time(sync.time() + length)

    def getMatrix(self, sync_part, part_index=-1):
        length = self.getLength(part_index)
        sync_part.time(sync_part.time() + length)
        matrix = np.eye(7)
        matrix[0, 6] = length
        return matrix


def make_action_lattice():
    def node(name, length):
        result = ActionTestNode(name)
        result.setLength(length)
        return result

    parent = node("parent", 0.6)
    parent.setnParts(2)
    entrance = node("entrance", 0.02)
    entrance.addChildNode(node("nested", 0.01), AccNode.EXIT)
    parent.addChildNode(entrance, AccNode.ENTRANCE)
    parent.addChildNode(node("before", 0.03), AccNode.BODY, 0, AccNode.BEFORE)
    parent.addChildNode(node("after", 0.04), AccNode.BODY, 1, AccNode.AFTER)
    parent.addChildNode(node("exit", 0.05), AccNode.EXIT)
    return make_lattice([parent])


def capture_actions(events, beam_key):
    actions = AccActionsContainer()
    for place in (AccNode.ENTRANCE, AccNode.BODY, AccNode.EXIT):
        for order in (0, 1):
            def action(params, place=place, order=order):
                beam = params[beam_key]
                sync = beam.getSyncParticle() if beam_key == "bunch" else beam.sync_part
                x = beam.x(0) if beam_key == "bunch" else beam.centroid[0]
                node = params["node"]
                assert params["actions"] is actions
                assert params["custom"] == "retained"
                events.append((node.getName(), params["parentNode"].getName(),
                               node.getActivePartIndex(), place, order,
                               params["path_length"], sync.time(), x))
                return object()  # Ignored by the action container.
            actions.addAction(action, place)
    return actions


@pytest.mark.parametrize("static", [False, True])
@pytest.mark.parametrize("fit", [False, True])
def test_envelope_actions_match_bunch_tracking(static, fit):
    lattice = make_action_lattice()
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    bunch.addParticle(0., 0., 0., 0., 0., 0.)
    envelope = Envelope(sync_part=bunch.getSyncParticle()).copy()
    expected = []
    bunch_actions = capture_actions(expected, "bunch")
    lattice.trackBunch(bunch, paramsDict={"custom": "retained", "path_length": 2.0},
                       actionContainer=bunch_actions)
    for _ in range(2):  # Cached traversal must replay every action.
        actual = []
        actions = capture_actions(actual, "envelope")
        original_actions = [list(actions.getActions(p)) for p in (0, 1, 2)]
        params = {"custom": "retained", "path_length": 2.0}
        tracked = envelope.copy()
        assert lattice.trackEnvelope(
            tracked, paramsDict=params, actionContainer=actions, static=static, fit=fit,
        ) is None
        assert len(actual) == len(expected)
        for observed, reference in zip(actual, expected):
            assert observed[:5] == reference[:5]
            assert np.allclose(observed[5:], reference[5:], atol=1e-10)
        assert params["envelope"] is tracked
        assert params["lattice"] is lattice
        assert params["path_length"] == pytest.approx(expected[-1][5])
        assert [actions.getActions(p) for p in (0, 1, 2)] == original_actions


@pytest.mark.parametrize("static", [False, True])
def test_envelope_empty_range_actions(static):
    lattice = make_action_lattice()
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    envelope = Envelope(sync_part=bunch.getSyncParticle())
    events = []
    actions = capture_actions(events, "envelope")
    params = {"custom": "retained", "path_length": 3.0}
    assert lattice.trackEnvelope(envelope, index_start=1, paramsDict=params,
                                  actionContainer=actions, static=static) is None
    assert events == []
    assert params["path_length"] == 3.0


def test_envelope_action_exception_preserves_container():
    lattice = make_action_lattice()
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    actions = AccActionsContainer()

    def fail(params):
        raise RuntimeError("Stop tracking")

    actions.addAction(fail, AccNode.BODY)
    with pytest.raises(RuntimeError, match="Stop tracking"):
        lattice.trackEnvelope(Envelope(sync_part=bunch.getSyncParticle()), actionContainer=actions)
    assert actions.getActions(AccNode.BODY) == [fail]


@pytest.mark.parametrize("fit", [False, True])
def test_envelope_body_actions_can_update_beam(fit):
    lattice = make_lattice([ActionTestNode("node")])
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    envelope = Envelope(sync_part=bunch.getSyncParticle())
    actions = AccActionsContainer()

    def change_beam(params):
        beam = params["envelope"]
        beam.centroid[0] += 0.2
        beam.sync_part.kinEnergy(beam.kin_energy + 0.1)

    actions.addAction(change_beam, AccNode.BODY)
    lattice.trackEnvelope(envelope, actionContainer=actions, fit=fit)
    assert envelope.centroid[0] == pytest.approx(0.2)
    assert envelope.kin_energy == pytest.approx(1.4)


@pytest.mark.parametrize("static", [False, True])
def test_envelope_uses_node_traversal_overrides(static):
    class CustomTraversalNode(ActionTestNode):
        def __init__(self):
            super().__init__("custom")
            self.visits = 0

        def trackActions(self, actionsContainer, paramsDict):
            self.visits += 1
            super().trackActions(actionsContainer, paramsDict)

    node = CustomTraversalNode()
    lattice = make_lattice([node])
    node.visits = 0  # Lattice initialization also uses the node traversal.
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    envelope = Envelope(sync_part=bunch.getSyncParticle())
    actions = monitor_actions(EnvelopeMonitor())
    for _ in range(2):
        lattice.trackEnvelope(envelope.copy(), actionContainer=actions, static=static)
    assert node.visits == 2


def test_envelope_cached_matrices_and_combined_fast_path(monkeypatch):
    node = ActionTestNode("cached")
    node.setLength(0.1)
    lattice = make_lattice([node])
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    envelope = Envelope(sync_part=bunch.getSyncParticle())
    monitor = EnvelopeMonitor()
    actions = monitor_actions(monitor)
    lattice.trackEnvelope(envelope.copy(), static=True, fit=False, actionContainer=actions)

    def unexpected_matrix(*args, **kwargs):
        pytest.fail("Cached tracking should reuse the existing node matrix")

    monkeypatch.setattr(node, "getMatrix", unexpected_matrix)
    lattice.trackEnvelope(envelope.copy(), static=True, fit=False, actionContainer=actions)
    assert len(monitor.history["s"]) == 4

    # Reuse the node cache to build the combined matrix, then use its fast path.
    expected = envelope.copy()
    lattice.trackEnvelope(expected, static=True, fit=False)
    tracked = envelope.copy()
    params = {"path_length": 2.0}
    lattice.trackEnvelope(tracked, static=True, fit=False, paramsDict=params)
    assert np.allclose(tracked.centroid, expected.centroid)
    assert np.allclose(tracked.cov_matrix, expected.cov_matrix)
    assert tracked.sync_part.time() == expected.sync_part.time()
    assert params["path_length"] == pytest.approx(2.1)
    assert params["node"] is node


@pytest.mark.parametrize("static", [False, True])
def test_static_time_advances_across_reused_matrices(static):
    lattice = make_lattice([ActionTestNode("periodic")])
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    envelope = Envelope(sync_part=bunch.getSyncParticle())
    first = envelope.copy()
    lattice.trackEnvelope(first, static=static, fit=False)
    initial_time = first.sync_part.time()
    lattice.trackEnvelope(first, static=static, fit=False)

    reference = envelope.copy()
    lattice.trackEnvelope(reference, fit=False)
    lattice.trackEnvelope(reference, fit=False)
    assert first.sync_part.time() == pytest.approx(2 * initial_time)
    assert first.sync_part.time() == pytest.approx(reference.sync_part.time())
    assert np.allclose(first.centroid, reference.centroid)
    assert np.allclose(first.cov_matrix, reference.cov_matrix)
