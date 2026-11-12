import numpy as np
import pytest

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.linac import MatrixRfGap
from orbit.bunch_utils import collect_bunch
from orbit.envelope import Envelope
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
    particles = np.random.default_rng(42).multivariate_normal(np.zeros(6), cov_matrix, size=nparts)
    for x in particles:
        bunch.addParticle(*x)

    # Covariance matrix calculated from particles will be slightly different.
    cov_matrix = calc_bunch_cov(bunch)

    envelope = Envelope(bunch=bunch).copy()
    data["bunch"]["cov"]["in"] = cov_scale * calc_bunch_cov(bunch)
    lattice.trackBunch(bunch)
    data["bunch"]["cov"]["out"] = cov_scale * calc_bunch_cov(bunch)

    # Track envelope
    analysis = BunchTwissAnalysis()
    analysis.analyzeBunch(bunch)
    bunch_centroid = np.array([analysis.getAverage(i) for i in range(6)])

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

    output_scale = np.sqrt(np.diag(data["bunch"]["cov"]["out"]))
    covariance_error = data["env"]["cov"]["out"] - data["bunch"]["cov"]["out"]
    assert np.max(np.abs(covariance_error) / np.outer(output_scale, output_scale)) < 0.01
    assert np.max(np.abs(envelope.centroid - bunch_centroid) / (output_scale / 1000.0)) < 0.01


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
    matrix_analytic = lattice.trackEnvelope(envelope.copy(), fit=False, return_matrix=True)
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
    matrix_actual = lattice.trackEnvelope(envelope.copy(), return_matrix=True, fit=True)
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
    # A cold uniform sphere obeys R'' = C/R**2, C = N*r0/(beta*gamma)**2.
    # Integrating this equation gives the distance at which its radius doubles.
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(0.01)
    radius = 0.01
    gamma = bunch.getSyncParticle().gamma()
    covariance = np.diag([radius**2 / 5, 0, radius**2 / 5, 0, radius**2 / (5 * gamma**2), 0])
    envelope = Envelope(sync_part=bunch.getSyncParticle(), cov_matrix=covariance, intensity=1e8)
    coefficient = envelope.intensity * envelope.classical_radius / (envelope.beta * gamma)**2
    distance = np.sqrt(radius**3 / (2 * coefficient)) * (np.sqrt(2) + np.arcsinh(1))
    lattice = make_lattice([DriftTEAPOT(length=distance, nparts=1000)])
    lattice.trackEnvelope(envelope, sc="3d", fit=False, static=True)
    radii = np.sqrt(5 * np.diag(envelope.cov_matrix)[[0, 2, 4]]) * [1, 1, gamma]
    np.testing.assert_allclose(radii, 2 * radius, rtol=0.002)


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
@pytest.mark.parametrize("static", [False, True])
def test_get_total_matrix(charge: float, static: bool) -> None:
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

    matrix = lattice.trackEnvelope(envelope.copy(), return_matrix=True, sc="2d", static=static)
    envelope_out_b = envelope.copy()
    envelope_out_b.transform(matrix)
    assert np.all(np.isclose(envelope_out_a.cov_matrix, envelope_out_b.cov_matrix))


def make_test_envelope():
    bunch = Bunch()
    bunch.getSyncParticle().kinEnergy(1.3)
    bunch.getSyncParticle().time(10.0)
    return Envelope(sync_part=bunch.getSyncParticle(), cov_matrix=make_default_cov_matrix()).copy()


@pytest.mark.parametrize("fit", [False, True])
def test_static_reuse_options_and_reset(monkeypatch, fit):
    from orbit.envelope import EnvelopeTracker

    lattice = make_lattice([DriftTEAPOT(length=0.2), DriftTEAPOT(length=0.7)])
    envelope = make_test_envelope()
    calls = []
    original = EnvelopeTracker.precompute

    def precompute(self, envelope, params_dict=None):
        calls.append((self.index_start, self.index_stop, self.fit, self.fallback))
        return original(self, envelope, params_dict=params_dict)

    monkeypatch.setattr(EnvelopeTracker, "precompute", precompute)
    expected = envelope.copy()
    for _ in range(3):
        lattice.trackEnvelope(expected, fit=fit)
    for _ in range(3):
        lattice.trackEnvelope(envelope, static=True, fit=fit)
    assert len(calls) == 1
    assert envelope.sync_part.time() == pytest.approx(expected.sync_part.time(), rel=0, abs=1e-14)
    np.testing.assert_allclose(envelope.cov_matrix, expected.cov_matrix)
    for options in [dict(index_stop=0, fit=fit), dict(index_stop=0, fit=not fit),
                    dict(index_stop=0, fit=not fit, fallback=False)]:
        actual = make_test_envelope()
        reference = actual.copy()
        reference_lattice = make_lattice([DriftTEAPOT(length=0.2), DriftTEAPOT(length=0.7)])
        reference_lattice.trackEnvelope(reference, **options)
        lattice.trackEnvelope(actual, static=True, **options)
        np.testing.assert_allclose(actual.cov_matrix, reference.cov_matrix)
    assert len(calls) == 4
    lattice.initialize()
    lattice.trackEnvelope(envelope, static=True)
    assert len(calls) == 5
    lattice.trackEnvelope(envelope, static=False)
    lattice.trackEnvelope(envelope, static=True)
    assert len(calls) == 6


@pytest.mark.parametrize("static", [False, True])
@pytest.mark.parametrize("start,stop,length", [(0, 0, .1), (1, 1, .2), (1, 2, .5),
    (0, -1, .6), (0, -2, .6), (-1, None, .6), (3, None, 0), (2, 0, 0)])
def test_sublattice_results(static, start, stop, length):
    from orbit.utils.matrix import get_matrix_drift

    lattice = make_lattice([DriftTEAPOT(length=x) for x in [.1, .2, .3]])
    envelope = make_test_envelope()
    reference = envelope.copy()
    expected_matrix = get_matrix_drift(reference.sync_part, length)
    matrix = lattice.trackEnvelope(envelope, index_start=start, index_stop=stop,
                                   static=static, fit=False, return_matrix=True)
    np.testing.assert_allclose(matrix, expected_matrix, atol=1e-14)
    assert envelope.sync_part.time() == pytest.approx(reference.sync_part.time(), rel=0, abs=1e-14)


@pytest.mark.parametrize("static", [False, True])
@pytest.mark.parametrize("fit", [False, True])
def test_actions_match_bunch_traversal(static, fit):
    parent = DriftTEAPOT(name="parent", length=.4, nparts=2)
    for place, part, location in [(AccNode.ENTRANCE, 0, AccNode.AFTER),
        (AccNode.BODY, 0, AccNode.BEFORE), (AccNode.BODY, 1, AccNode.AFTER),
        (AccNode.EXIT, 0, AccNode.AFTER)]:
        parent.addChildNode(DriftTEAPOT(name=f"child-{place}-{part}", length=.1),
                            place, part_index=part, place_in_part=location)
    lattice = make_lattice([parent])
    envelope = make_test_envelope()
    envelope.centroid[1] = .001
    bunch = Bunch()
    bunch.mass(envelope.mass)
    bunch.charge(envelope.charge)
    bunch.getSyncParticle().kinEnergy(envelope.kin_energy)
    bunch.getSyncParticle().time(10)
    bunch.addParticle(*envelope.centroid)

    def observe(beam_key):
        actions = AccActionsContainer()
        events = []
        for place in [AccNode.ENTRANCE, AccNode.BODY, AccNode.EXIT]:
            def action(params_dict, place=place):
                assert params_dict["actions"] is actions
                assert params_dict["custom"] == "retained"
                beam = params_dict[beam_key]
                sync_part = beam.sync_part if beam_key == "envelope" else beam.getSyncParticle()
                x = beam.centroid[0] if beam_key == "envelope" else beam.x(0)
                events.append((params_dict["node"].getName(), params_dict["parentNode"].getName(),
                    place, params_dict["node"].getActivePartIndex(), params_dict["path_length"],
                    sync_part.time(), x))
            actions.addAction(action, place)
        return actions, events

    bunch_actions, expected = observe("bunch")
    lattice.trackBunch(bunch, paramsDict={"custom": "retained"}, actionContainer=bunch_actions)
    actions, actual = observe("envelope")
    original = [list(actions.getActions(place)) for place in range(3)]
    lattice.trackEnvelope(envelope, params_dict={"custom": "retained"}, action_container=actions,
                          static=static, fit=fit)
    assert len(actual) == len(expected)
    for observed, reference in zip(actual, expected):
        assert observed[:4] == reference[:4]
        np.testing.assert_allclose(observed[4:], reference[4:], rtol=0, atol=1e-12)
    assert [actions.getActions(place) for place in range(3)] == original


@pytest.mark.parametrize("failure", ["action", "matrix"])
def test_actions_removed_after_failure(failure):
    class FailingDrift(DriftTEAPOT):
        def getMatrix(self, sync_part, part_index=-1):
            raise RuntimeError("matrix failure")

    node = FailingDrift(length=.1) if failure == "matrix" else DriftTEAPOT(length=.1)
    lattice = make_lattice([node])
    actions = AccActionsContainer()
    def fail(params_dict):
        raise RuntimeError("action failure")
    if failure == "action":
        actions.addAction(fail, AccNode.EXIT)
    original = [list(actions.getActions(place)) for place in range(3)]
    with pytest.raises(RuntimeError, match=f"{failure} failure"):
        lattice.trackEnvelope(make_test_envelope(), fit=False, action_container=actions)
    assert [actions.getActions(place) for place in range(3)] == original
    if failure == "action":
        actions.removeAction(fail, AccNode.EXIT)
    healthy = make_lattice([DriftTEAPOT(length=.1)])
    healthy.trackEnvelope(make_test_envelope(), action_container=actions)
    assert all(not actions.getActions(place) for place in range(3))


@pytest.mark.parametrize("static", [False, True])
@pytest.mark.parametrize("fallback", [False, True])
def test_missing_matrix_fallback(static, fallback):
    class UnsupportedDrift(DriftTEAPOT):
        def getMatrix(self, sync_part, part_index=-1):
            raise NotImplementedError("missing matrix")
    lattice = make_lattice([UnsupportedDrift(length=.3)])
    envelope = make_test_envelope()
    if not fallback:
        with pytest.raises(NotImplementedError, match="missing matrix"):
            lattice.trackEnvelope(envelope, fit=False, fallback=False, static=static)
    else:
        expected = envelope.copy()
        lattice.trackEnvelope(expected, fit=True)
        lattice.trackEnvelope(envelope, fit=False, fallback=True, static=static)
        np.testing.assert_allclose(envelope.cov_matrix, expected.cov_matrix)


@pytest.mark.parametrize("result", [None, ValueError("invalid matrix")])
def test_fallback_respects_none_and_other_errors(result):
    class CustomNode(AccNode):
        def getMatrix(self, sync_part, part_index=-1):
            if isinstance(result, Exception):
                raise result
            return result
        def track(self, params_dict):
            pytest.fail("Unexpected fitting")
    lattice = make_lattice([CustomNode()])
    envelope = make_test_envelope()
    if isinstance(result, Exception):
        with pytest.raises(ValueError, match="invalid matrix"):
            lattice.trackEnvelope(envelope, fit=False)
    else:
        expected = envelope.copy()
        lattice.trackEnvelope(envelope, fit=False)
        np.testing.assert_array_equal(envelope.cov_matrix, expected.cov_matrix)


@pytest.mark.parametrize("fit", [False, True])
def test_reference_updates_and_preparation_independence(fit):
    from orbit.envelope import EnvelopeTracker
    class AcceleratingNode(AccNode):
        def update(self, sync_part):
            sync_part.kinEnergy(sync_part.kinEnergy() + .1)
            sync_part.time(sync_part.time() + .25)
            sync_part.rVector(tuple(x + .1 for x in sync_part.rVector()))
            sync_part.nxVector((0, 1, 0))
        def getMatrix(self, sync_part, part_index=-1):
            self.update(sync_part)
            return np.eye(7)
        def track(self, params_dict):
            self.update(params_dict["bunch"].getSyncParticle())
    lattice = make_lattice([AcceleratingNode(), AcceleratingNode()])
    envelope = make_test_envelope()
    tracker = EnvelopeTracker(lattice, fit=fit)
    params_dict = {"path_length": 5}
    tracker.precompute(envelope, params_dict=params_dict)
    assert params_dict == {"path_length": 5}
    assert envelope.sync_part.time() == 10
    assert envelope.kin_energy == pytest.approx(1.3)
    np.testing.assert_array_equal(envelope.sync_part.rVector(), (0, 0, 0))
    # Dynamic transport must recalculate at the updated energy on each turn.
    lattice.trackEnvelope(envelope, fit=fit)
    lattice.trackEnvelope(envelope, fit=fit)
    assert envelope.kin_energy == pytest.approx(1.7)
    assert envelope.sync_part.time() == pytest.approx(11)
    np.testing.assert_allclose(envelope.sync_part.rVector(), (.4, .4, .4))
    np.testing.assert_allclose(envelope.sync_part.nxVector(), (0, 1, 0))


def test_envelope_copy_owns_reference_state():
    import gc
    original = make_test_envelope()
    original.sync_part.rVector((1, 2, 3))
    original.sync_part.nxVector((0, 1, 0))
    original.centroid[:] = np.arange(6)
    original.intensity = 123
    original.rms_bunch_length_factor = 4
    copied = original.copy()
    for attr in ["mass", "charge", "kinEnergy", "time", "rVector", "pVector", "nxVector"]:
        np.testing.assert_allclose(getattr(copied.sync_part, attr)(), getattr(original.sync_part, attr)())
    assert copied.intensity == original.intensity
    assert copied.classical_radius == original.classical_radius
    assert copied.rms_bunch_length_factor == 4
    copied.centroid[0] = 99
    copied.cov_matrix[0, 0] = 99
    copied.sync_part.time(20)
    copied.sync_part.kinEnergy(.5)
    assert original.centroid[0] == 0
    assert original.cov_matrix[0, 0] != 99
    assert original.sync_part.time() == 10
    assert original.kin_energy == pytest.approx(1.3)
    del original
    gc.collect()
    assert copied.copy().kin_energy == pytest.approx(.5)


@pytest.mark.parametrize("kind", ["teapot", "linac"])
@pytest.mark.parametrize("fringe", [False, True])
@pytest.mark.parametrize("static", [False, True])
@pytest.mark.parametrize("ea1,ea2", [(0.1, 0), (0, -.05), (.1, -.05)])
def test_edge_bends_with_fringe(kind, fringe, static, ea1, ea2):
    if kind == "teapot":
        node = BendTEAPOT(length=1, theta=.2, ea1=ea1, ea2=ea2, nparts=5)
    else:
        node = Bend()
        node.setLength(1)
        node.setnParts(5)
        for name, value in [("theta", .2), ("ea1", ea1), ("ea2", ea2)]:
            node.setParam(name, value)
    lattice = make_lattice([node])
    if hasattr(node, "setUsageFringeFieldIN"):
        node.setUsageFringeFieldIN(fringe)
        node.setUsageFringeFieldOUT(fringe)
    envelope = make_test_envelope()
    bunch = Bunch()
    bunch.mass(envelope.mass)
    bunch.charge(envelope.charge)
    bunch.getSyncParticle().kinEnergy(envelope.kin_energy)
    expected = fit_transfer_matrix(lattice, bunch)
    actual = lattice.trackEnvelope(envelope, fit=False, static=static, return_matrix=True)
    np.testing.assert_allclose(actual, expected, atol=1e-10)
    with pytest.raises(NotImplementedError, match="edge angles"):
        node.getMatrix(envelope.sync_part, part_index=0)
    for fringe_node in (node.getNodeFringeFieldIN(), node.getNodeFringeFieldOUT()):
        with pytest.raises(NotImplementedError, match="edge fringes"):
            fringe_node.getMatrix(envelope.sync_part)
    with pytest.raises(NotImplementedError):
        lattice.trackEnvelope(envelope.copy(), fit=False, fallback=False)
