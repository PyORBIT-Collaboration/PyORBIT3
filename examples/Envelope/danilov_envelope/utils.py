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
    fill_frac: float = 0.5,
    start: str = "drift",
    fringe: bool = False,
    max_part_length: float = 0.1,
    verbose: bool = False,
) -> AccLattice:
    """Create FODO lattice with specified phase advances.

    Args:
        phase_adv_x: Phase advance in x plane [rad].
        phase_adv_y: Phase advance in y plane [rad].
        length: Length of the lattice [m].
        mass: Particle mass [GeV/c^2]
        kin_energy: Synchronous particle kinetic energy [GeV].
        fill_frac : Fraction of the lattice occupied by quadrupoles.
        start : str
            "drift":
                O-F-F-O-O-D-D-O. 
            "quad":
                F-O-O-D-D-O-O-F
        fringe: Toggles fringe fields before/after quads.
        max_part_length: Maximum node part length [m].
        verbose: Print optimization status.

    Returns:
        TEAPOT_Lattice
    """

    def _make_lattice(k1: float, k2: float) -> AccLattice:
        """Create FODO lattice with specified focusing strengths.

        k1 and k2 are the focusing strengths of the
        focusing (1st) and defocusing (2nd) quads, respectively.
        """
        length_quad = length * fill_frac / 2.0
        length_drift = length * (1.0 - fill_frac) / 2.0

        if start == "quad":
            drift_nodes = [
                DriftTEAPOT("drift1"),
                DriftTEAPOT("drift2"),
            ]
            quad_nodes = [
                QuadTEAPOT("qf1"),
                QuadTEAPOT("qd"),
                QuadTEAPOT("qf2"),
            ]

            drift_nodes[0].setLength(length_drift)
            drift_nodes[1].setLength(length_drift)

            quad_nodes[0].setLength(length_quad * 0.5)
            quad_nodes[1].setLength(length_quad)
            quad_nodes[2].setLength(length_quad * 0.5)

            quad_nodes[0].setParam("kq", +k1)
            quad_nodes[1].setParam("kq", -k2)
            quad_nodes[2].setParam("kq", +k1)

            lattice = TEAPOT_Lattice()
            lattice.addNode(quad_nodes[0])
            lattice.addNode(drift_nodes[0])
            lattice.addNode(quad_nodes[1])
            lattice.addNode(drift_nodes[1])
            lattice.addNode(quad_nodes[2])
            lattice.initialize()

        elif start == "drift":

            drift_nodes = [
                DriftTEAPOT("drift1"),
                DriftTEAPOT("drift2"),
                DriftTEAPOT("drift3"),
            ]
            quad_nodes = [
                QuadTEAPOT("qf"),
                QuadTEAPOT("qd"),
            ]

            drift_nodes[0].setLength(length_drift * 0.5)
            drift_nodes[1].setLength(length_drift)
            drift_nodes[2].setLength(length_drift * 0.5)

            quad_nodes[0].setLength(length_quad)
            quad_nodes[1].setLength(length_quad)

            quad_nodes[0].setParam("kq", +k1)
            quad_nodes[1].setParam("kq", -k2)

            lattice = TEAPOT_Lattice()
            lattice.addNode(drift_nodes[0])
            lattice.addNode(quad_nodes[0])
            lattice.addNode(drift_nodes[1])
            lattice.addNode(quad_nodes[1])
            lattice.addNode(drift_nodes[2])
            lattice.initialize()

            # Toggle fringe fields
            for node in lattice.getNodes():
                node.setUsageFringeFieldIN(fringe)
                node.setUsageFringeFieldOUT(fringe)

            lattice.initialize()
            return lattice

    def loss_function(parameters: np.ndarray) -> float:
        lattice = _make_lattice(*parameters)
        phase_adv_calc = get_phase_adv(lattice, mass=mass, kin_energy=kin_energy)
        phase_adv_targ = np.array([phase_adv_x, phase_adv_y])
        return np.abs(phase_adv_calc - phase_adv_targ)

    parameters = np.array([0.5, 0.5])  # ~ 80 deg phase advance
    result = scipy.optimize.least_squares(loss_function, parameters, verbose=verbose)
    lattice = _make_lattice(*result.x)

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
