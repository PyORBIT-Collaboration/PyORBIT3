import copy
import numpy as np
import scipy.optimize

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist2D
from orbit.bunch_generators import KVDist2D
from orbit.bunch_generators import WaterBagDist2D
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT


def split_node(node: AccNode, max_part_length: float = None) -> AccNode:
    if max_part_length:
        if node.getLength() > max_part_length:
            node.setnParts(1 + int(node.getLength() / max_part_length))
    return node


def get_phase_adv(lattice: AccLattice, mass: float, kin_energy: float) -> np.ndarray:
    bunch = Bunch()
    bunch.mass(mass)
    bunch.getSyncParticle().kinEnergy(kin_energy)
    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
    lattice_params = matrix_lattice.getRingParametersDict()
    phase_adv = [
        lattice_params["fractional tune x"],
        lattice_params["fractional tune y"],
    ]
    phase_adv = np.array(phase_adv)
    phase_adv = phase_adv * 2.0 * np.pi
    return phase_adv


def make_fodo_lattice(
    phase_adv_x: float,
    phase_adv_y: float,
    length: float,
    mass: float,
    kin_energy: float,
    fill_factor: float = 0.5,
    start: str = "drift",
    fringe: bool = False,
    max_part_length: float = 0.1,
    verbose: bool = False,
) -> AccLattice:
    """Create FODO lattice with specified phase advances.

    Parameters
    ----------
    phase_adv_x{y}: float
        The x{y} lattice phase advance [rad].
    length : float
        The length of the lattice [m].
    mass, kin_energy : float
        Mass [GeV/c^2] and kinetic energy [GeV] of synchronous particle.
    fill_fac : float
        The fraction of the lattice occupied by quadrupoles.
    fringe : bool
        Whether to include nonlinear fringe fields in the lattice.
    start : str
        If 'drift', the lattice will be O-F-O-O-D-O. If 'quad' the lattice will
        be (F/2)-O-O-D-O-O-(F/2).
    reverse : bool
        If True, reverse the lattice elements.

    Returns
    -------
    TEAPOT_Lattice
    """

    def _make_lattice(k1: float, k2: float) -> AccLattice:
        """Create FODO lattice with specified focusing strengths.

        k1 and k2 are the focusing strengths of the
        focusing (1st) and defocusing (2nd) quads, respectively.
        """
        # Instantiate elements
        lattice = TEAPOT_Lattice()
        drift1 = DriftTEAPOT("drift1")
        drift2 = DriftTEAPOT("drift2")
        drift_half1 = DriftTEAPOT("drift_half1")
        drift_half2 = DriftTEAPOT("drift_half2")
        qf = QuadTEAPOT("qf")
        qd = QuadTEAPOT("qd")
        qf_half1 = QuadTEAPOT("qf_half1")
        qf_half2 = QuadTEAPOT("qf_half2")
        qd_half1 = QuadTEAPOT("qd_half1")
        qd_half2 = QuadTEAPOT("qd_half2")

        # Set lengths
        half_nodes = (drift_half1, drift_half2, qf_half1, qf_half2, qd_half1, qd_half2)
        full_nodes = (drift1, drift2, qf, qd)
        for node in half_nodes:
            node.setLength(length * fill_factor / 4.0)
        for node in full_nodes:
            node.setLength(length * fill_factor / 2.0)

        # Set quad focusing strengths
        for node in (qf, qf_half1, qf_half2):
            node.addParam("kq", +k1)
        for node in (qd, qd_half1, qd_half2):
            node.addParam("kq", -k2)

        # Create lattice
        if start == "drift":
            lattice.addNode(drift_half1)
            lattice.addNode(qf)
            lattice.addNode(drift2)
            lattice.addNode(qd)
            lattice.addNode(drift_half2)
        elif start == "quad":
            lattice.addNode(qf_half1)
            lattice.addNode(drift1)
            lattice.addNode(qd)
            lattice.addNode(drift2)
            lattice.addNode(qf_half2)

        # Toggle fringe fields
        for node in lattice.getNodes():
            node.setUsageFringeFieldIN(fringe)
            node.setUsageFringeFieldOUT(fringe)

        lattice.initialize()
        return lattice

    def function(k: np.ndarray) -> float:
        lattice = _make_lattice(k[0], k[1])
        phase_adv_calc = get_phase_adv(lattice, mass=mass, kin_energy=kin_energy)
        phase_adv_targ = np.array([phase_adv_x, phase_adv_y])
        return np.abs(phase_adv_calc - phase_adv_targ)

    k0 = np.array([0.5, 0.5])  # ~ 80 deg phase advance
    result = scipy.optimize.least_squares(function, k0, verbose=verbose)
    k1, k2 = result.x
    lattice = _make_lattice(k1, k2)

    for node in lattice.getNodes():
        node = split_node(node, max_part_length)

    if verbose:
        phase_adv_calc = get_phase_adv(lattice, mass=mass, kin_energy=kin_energy)
        phase_adv_targ = np.array([phase_adv_x, phase_adv_y])
        phase_adv_calc *= 180.0 / np.pi
        phase_adv_targ *= 180.0 / np.pi
        print(f"phase_adv_x = {phase_adv_calc[0]} (target={phase_adv_targ[0]})")
        print(f"phase_adv_y = {phase_adv_calc[1]} (target={phase_adv_targ[1]})")

    return lattice


def get_bunch_coords(bunch: Bunch, axis: tuple[int, ...] = None) -> np.ndarray:
    if axis is None:
        axis = tuple(range(6))

    X = np.zeros((bunch.getSize(), 6))
    for i in range(bunch.getSize()):
        X[i, 0] = bunch.x(i)
        X[i, 1] = bunch.xp(i)
        X[i, 2] = bunch.y(i)
        X[i, 3] = bunch.yp(i)
        X[i, 4] = bunch.z(i)
        X[i, 5] = bunch.dE(i)
    return X[:, axis]


def get_bunch_cov(bunch: Bunch) -> np.ndarray:
    calc = BunchTwissAnalysis()
    calc.computeBunchMoments(bunch, 2, 0, 0)

    cov_matrix = np.zeros((6, 6))
    for i in range(6):
        for j in range(i + 1):
            cov_matrix[i, j] = calc.getCorrelation(j, i)
            cov_matrix[j, i] = cov_matrix[i, j]
    return cov_matrix


class BunchMonitor:
    def __init__(self, verbose: int = 1) -> None:
        self.distance = 0.0
        self._pos_old = 0.0
        self._pos_new = 0.0
        self.verbose = verbose

        self.history = {}
        for key in [
            "s",
            "xrms",
            "yrms",
            "epsx",
            "epsy",
            "cov_00",
            "cov_01",
            "cov_02",
            "cov_03",
            "cov_11",
            "cov_12",
            "cov_13",
            "cov_22",
            "cov_23",
            "cov_33",
            "rxy",
        ]:
            self.history[key] = []

    def package_history(self) -> None:
        history = copy.deepcopy(self.history)
        for key in history:
            history[key] = np.array(history[key])
        history["s"] -= history["s"][0]
        return history

    def __call__(self, params_dict: dict) -> None:
        bunch = params_dict["bunch"]
        node = params_dict["node"]

        self._pos_new = params_dict["path_length"]
        if self._pos_old > self._pos_new:
            self._pos_old = 0.0
        self.distance += self._pos_new - self._pos_old
        self._pos_old = self._pos_new

        cov_matrix = get_bunch_cov(bunch)

        for i in range(4):
            for j in range(i, 4):
                key = f"cov_{i}{j}"
                self.history[key].append(cov_matrix[i, j])

        self.history["s"].append(self.distance)
        self.history["xrms"].append(np.sqrt(cov_matrix[0, 0]))
        self.history["yrms"].append(np.sqrt(cov_matrix[2, 2]))
        self.history["epsx"].append(np.sqrt(np.linalg.det(cov_matrix[0:2, 0:2])))
        self.history["epsy"].append(np.sqrt(np.linalg.det(cov_matrix[2:4, 2:4])))
        self.history["rxy"].append(
            self.history["cov_02"][-1]
            / np.sqrt(self.history["cov_00"][-1] * self.history["cov_22"][-1])
        )

        if self.verbose:
            message = ""
            message += "s={:0.3f} ".format(self.history["s"][-1])
            message += "xrms={:0.3f} ".format(self.history["xrms"][-1] * 1000.0)
            message += "yrms={:0.3f} ".format(self.history["yrms"][-1] * 1000.0)
            print(message)
