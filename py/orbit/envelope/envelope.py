import math
import warnings

import numpy as np
import scipy.constants
import scipy.special

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle
from orbit.lattice import AccNode
from orbit.lattice import AccLattice

from orbit.teapot import ApertureTEAPOT
from orbit.teapot import DriftTEAPOT
from orbit.teapot import BendTEAPOT
from orbit.teapot import KickTEAPOT
from orbit.teapot import MonitorTEAPOT
from orbit.teapot import MultipoleTEAPOT
from orbit.teapot import NodeTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import SolenoidTEAPOT
from orbit.teapot import FringeFieldTEAPOT
from orbit.teapot import BunchWrapTEAPOT
from orbit.teapot import TiltTEAPOT
from orbit.teapot import ContinuousLinearFocusingTEAPOT
from orbit.teapot import TurnCounterTEAPOT

from orbit.py_linac.lattice import MarkerLinacNode as MarkerLINAC
from orbit.py_linac.lattice import Drift as DriftLINAC
from orbit.py_linac.lattice import Quad as QuadLINAC
from orbit.py_linac.lattice import Bend as BendLINAC
from orbit.py_linac.lattice import DCorrectorH as DCorrectorHLINAC
from orbit.py_linac.lattice import DCorrectorV as DCorrectorVLINAC
from orbit.py_linac.lattice import Solenoid as SolenoidLINAC
from orbit.py_linac.lattice import TiltElement as TiltLINAC
from orbit.py_linac.lattice import FringeField as FringeFieldLINAC
from orbit.py_linac.lattice import BaseRF_Gap as BaseRF_Gap
from orbit.py_linac.lattice import LinacApertureNode as ApertureLINAC

from .matrix import get_dp_p_coeff
from .matrix import get_zp_coeff
from .matrix import convert_matrix_dp_p_to_dE
from .matrix import convert_matrix_zp_to_dE
from .matrix import track_sync_part_tilt
from .matrix import track_sync_part_kick
from .matrix import track_sync_part_drift
from .matrix import track_sync_part_quad
from .matrix import track_sync_part_bend
from .matrix import track_sync_part_solenoid
from .matrix import track_sync_part_rf_gap
from .matrix import track_sync_part_cf
from .utils import gen_dist
from .utils import proj_cov_matrix


ENTRANCE = AccNode.ENTRANCE
BODY = AccNode.BODY
EXIT = AccNode.EXIT

BEFORE = AccNode.BEFORE
AFTER = AccNode.AFTER

IGNORE_NODE_TYPES = [
    NodeTEAPOT,
    MonitorTEAPOT,
    FringeFieldTEAPOT,
    ApertureTEAPOT,
    BunchWrapTEAPOT,
    TurnCounterTEAPOT,
    MarkerLINAC,
    FringeFieldLINAC,
]


def build_diag_matrix_from_xyz_eig(eigenvectors: np.ndarray) -> np.ndarray:
    A = np.eye(7)
    for i in range(eigenvectors.shape[0]):
        for j in range(eigenvectors.shape[1]):
            row = i * 2
            col = j * 2
            A[row, col] = A[row + 1, col + 1] = eigenvectors[i, j]
    return A


def track_sync_part(
    node: AccNode,
    sync_part: SyncParticle,
    charge: float,
    index: int = -1,
) -> np.ndarray | None:
    """Calculate transfer matrix and update synchronous particle.

    This function maps various accelerator nodes to 7 x 7 transfer matrices
    for envelope tracking. For non-accelerating, finite-length nodes, the
    synchronous particle time is updated as in a drift. Accelerating nodes
    such as RF gaps will update the synchronous particle energy.

    Args:
        node: The accelerator node.
        sync_part: Synchronous particle.
        charge: Particle charge. (The charge is currently an attribute of the
            bunch, not the synchronous particle.)
        index: Node part index. An index of -1 will return the transfer matrix
            for the entire node.
    Returns:
        7 x 7 transfer matrix or None. If None, the node can be ignored during
        envelope tracking.
    """
    node_type = type(node)
    if node_type in IGNORE_NODE_TYPES:
        return None

    length = node.getLength(index)
    nparts = node.getnParts()

    if node_type is DriftTEAPOT:
        if length <= 0:
            return None
        return track_sync_part_drift(sync_part=sync_part, length=length)

    elif node_type is SolenoidTEAPOT:
        if length <= 0:
            return None
        B = node.getParam("B")
        if node.waveform:
            B *= self.waveform.getStrength()
        return track_sync_part_solenoid(sync_part=sync_part, length=length, B=B, charge=charge)

    elif node_type is MultipoleTEAPOT:
        if length <= 0:
            return None
        if np.all(np.abs(node.getParam("kls")) == 0):
            return track_sync_part_drift(sync_part=sync_part, length=length)

    elif node_type is QuadTEAPOT:
        if length <= 0:
            return None
        kq = node.getParam("kq")
        if node.waveform:
            kq *= node.waveform.getStrength()
        return track_sync_part_quad(sync_part=sync_part, length=length, kq=kq, charge=charge)

    elif node_type is BendTEAPOT:
        if length <= 0:
            return None
        theta = node.getParam("theta") / (nparts - 1)
        if index == 0 or index == nparts - 1:
            theta *= 0.5
        return track_sync_part_bend(sync_part=sync_part, length=length, theta=theta, charge=charge)

    elif node_type is KickTEAPOT:
        scale = 1.0
        if node.waveform is not None:
            scale = node.waveform.getStrength()

        scale /= (nparts - 1)
        kx = scale * node.getParam("kx")
        ky = scale * node.getParam("ky")
        kE = node.getParam("dE")

        if abs(kx) > 0 or abs(ky) > 0 or abs(kE) > 0:
            return np.matmul(
                track_sync_part_kick(sync_part=sync_part, kx=kx, ky=ky, kE=kE),
                track_sync_part_drift(sync_part=sync_part, length=length),
            )
        else:
            return track_sync_part_drift(sync_part=sync_part, length=length)

    elif node_type is TiltTEAPOT:
        angle = node.getTiltAngle()
        if angle == 0:
            return None
        return track_sync_part_tilt(sync_part=sync_part, angle=angle)

    elif node_type is ContinuousLinearFocusingTEAPOT:
        if length <= 0:
            return None
        kq = node.getParam("kq")
        if node.waveform:
            kq *= node.waveform.getStrength()
        return track_sync_part_cf(sync_part=sync_part, length=length, kq=kq)

    elif node_type is DriftLINAC:
        if length <= 0:
            return None
        return track_sync_part_drift(sync_part=sync_part, length=length)

    elif node_type is QuadLINAC:
        if length <= 0:
            return None
        brho = 3.335640952 * sync_part.momentum() / charge
        kq = node.getParam("dB/dr") / brho
        return track_sync_part_quad(sync_part=sync_part, length=length, kq=kq, charge=charge)

    elif node_type is BendLINAC:
        if length <= 0:
            return None
        theta = node.getParam("theta") / (nparts - 1)
        if index == 0 or index == nparts - 1:
            theta *= 0.5
        return track_sync_part_bend(sync_part=sync_part, length=length, theta=theta, charge=charge)

    elif node_type is DCorrectorHLINAC:
        length = node.getParam("effLength") / nparts
        field = node.getParam("B")
        delta_xp = -field * charge * length * 0.299792 / sync_part.momentum()
        if delta_xp == 0:
            return None
        return track_sync_part_kick(sync_part=sync_part, kx=delta_xp, ky=0.0, kE=0.0)

    elif node_type is DCorrectorVLINAC:
        length = node.getParam("effLength") / nparts
        field = node.getParam("B")
        delta_yp = -field * charge * length * 0.299792 / sync_part.momentum()
        if delta_yp == 0:
            return None
        return track_sync_part_kick(sync_part=sync_part, kx=0.0, ky=delta_yp, kE=0.0)

    elif node_type is SolenoidLINAC:
        if length <= 0:
            return None
        B = node.getParam("B")
        return track_sync_part_solenoid(sync_part=sync_part, length=length, B=B, charge=charge)

    elif node_type is TiltLINAC:
        angle = node.getTiltAngle()
        if angle == 0:
            return None
        return track_sync_part_tilt(sync_part=sync_part, angle=angle)

    elif node_type is BaseRF_Gap:
        E0TL = node.getParam("E0TL")
        mode_phase = node.getParam("mode") * math.pi

        cavity = node.getRF_Cavity()
        frequency = cavity.getFrequency()
        phase = cavity.getPhase() + mode_phase
        amplitude = cavity.getAmp()

        arrival_time = sync_part.time()
        arrival_time_design = cavity.getDesignArrivalTime()

        if node.isFirstRFGap():
            if cavity.isDesignSetUp():
                phase = math.fmod(frequency * (arrival_time - arrival_time_design) * 2.0 * math.pi + phase, 2.0 * math.pi)
            else:
                orbitFinalize("Run `trackDesign` first to initialize cavity phases.")
        else:
            phase = math.fmod(frequency * (arrival_time - arrival_time_design) * 2.0 * math.pi + phase,2.0 * math.pi)

        node.setGapPhase(phase)

        if amplitude == 0.0:
            return None

        return track_sync_part_rf_gap(
            sync_part=sync_part,
            frequency=frequency,
            E0TL=(E0TL * amplitude),
            phase=phase,
            charge=charge,
        )

    raise NotImplementedError(str(node))


class Envelope:
    """Represents beam envelope/centroid.

    Attributes:
        bunch: Bunch containing synchronous particle and (optionally) test particles.
        cov_matrix: 6 x 6 covariance matrix
        centroid: 6 x 1 centroid vector.
        intensity: Total number of particles.
    """

    def __init__(
        self,
        bunch: Bunch,
        cov_matrix: np.ndarray = None,
        centroid: np.ndarray = None,
        intensity: float = 0.0,
    ) -> None:

        # Eventually allow:
        #   - setting covariance matrix from bunch particles
        #   - tracking bunch particles as test particles
        empty_bunch = Bunch()
        bunch.copyEmptyBunchTo(empty_bunch)

        self.bunch = empty_bunch
        self.sync_part = empty_bunch.getSyncParticle()

        self.classical_radius = self.bunch.classicalRadius()

        self.centroid = centroid
        if self.centroid is None:
            self.centroid = np.zeros(6)

        self.cov_matrix = cov_matrix
        if self.cov_matrix is None:
            self.cov_matrix = np.eye(6)

        self.intensity = 0.0
        self.set_intensity(intensity)

        self.rms_bunch_length_factor = np.sqrt(12.0)

    def set_intensity(self, intensity: float) -> None:
        self.intensity = intensity

    @property
    def sc_factor(self) -> float:
        return (
            2.0
            * self.intensity
            * self.classical_radius
            / (self.beta() ** 2 * self.gamma() ** 3)
        )

    def gamma(self) -> float:
        return self.sync_part.gamma()

    def beta(self) -> float:
        return self.sync_part.beta()

    def mass(self) -> float:
        return self.sync_part.mass()

    def charge(self) -> float:
        return self.bunch.charge()

    def rms(self, axis: int = None) -> float | np.ndarray:
        rms_arr = np.sqrt(np.diag(self.cov_matrix))
        return rms_arr[axis]

    def transform(self, matrix: np.ndarray) -> None:
        m = matrix[:-1, :-1]
        u = matrix[:-1, -1]
        self.cov_matrix = m @ self.cov_matrix @ m.T
        self.centroid = np.matmul(m, self.centroid) + u

    def sample(self, size: int, dist: str = "kv") -> np.ndarray:
        particles = gen_dist(size=size, cov_matrix=self.cov_matrix, name=dist)
        particles = particles + self.centroid
        return particles

    def sc_matrix_2d(self, length: float) -> np.ndarray:
        centroid = self.centroid
        cov_matrix = self.cov_matrix

        # Calculate transfer matrix in normalized (upright) frame.
        cov_xx = cov_matrix[0, 0]
        cov_yy = cov_matrix[2, 2]
        cov_xy = cov_matrix[0, 2]

        phi = -0.5 * np.arctan2(2 * cov_xy, cov_xx - cov_yy)
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)
        rx = 2.0 * np.sqrt(abs(cov_xx * cos_phi**2 + cov_yy * sin_phi**2 - 2.0 * cov_xy * sin_phi * cos_phi))
        ry = 2.0 * np.sqrt(abs(cov_xx * sin_phi**2 + cov_yy * cos_phi**2 + 2.0 * cov_xy * sin_phi * cos_phi))

        bunch_length = self.rms_bunch_length_factor * np.sqrt(cov_matrix[4, 4])
        perveance = self.sc_factor / bunch_length
        kappa_factor = 2.0 * perveance / (rx + ry)

        M = np.identity(7)
        M[1, 0] = kappa_factor * length / rx
        M[3, 2] = kappa_factor * length / ry

        # Build matrix A to transform out of normalized frame.
        A = np.eye(7)
        A[0, 0] = A[1, 1] = +cos_phi
        A[0, 2] = A[1, 3] = +sin_phi
        A[2, 0] = A[3, 1] = -sin_phi
        A[2, 2] = A[3, 3] = +cos_phi

        A_inv = A.T

        # Build matrix T to shift to beam centroid.
        T = np.identity(7)
        T[0, -1] = centroid[0]
        T[2, -1] = centroid[2]

        T_inv = np.copy(T)
        T_inv[:-1, -1] = -T[:-1, -1]

        # Compute transfer matrix in lab frame.
        return T @ A @ M @ A_inv @ T_inv

    def sc_matrix_3d(self, length: float) -> np.ndarray:
        # Build Lorentz matrix: rest frame to lab frame.
        # x -> x
        # y -> y
        # z -> z / gamma
        # x' = dx/ds -> x' * gamma
        # y' = dy/ds -> y' * gamma
        # z' = dz/ds -> z'
        gamma = self.gamma()
        gamma_inv = 1.0 / gamma

        L = np.identity(7)
        L[1, 1] = gamma
        L[3, 3] = gamma
        L[4, 4] = gamma_inv
        L_inv = np.diag(1.0 / np.diag(L))

        # Get centroid in rest frame.
        centroid = np.matmul(L_inv[:-1, :-1], self.centroid)

        # Get covariance matrix in rest frame.
        cov_matrix = L_inv[:-1, :-1] @ self.cov_matrix @ L_inv[:-1, :-1].T

        # Project covariance matrix onto x-y-z plane.
        cov_matrix_proj = proj_cov_matrix(cov_matrix, axis=(0, 2, 4))

        # Compute eigenvalues and eigenvectors of x-y-z covariance matrix.
        cov_eig_res = np.linalg.eigh(cov_matrix_proj)
        cov_eig_vals = cov_eig_res.eigenvalues
        cov_eig_vecs = cov_eig_res.eigenvectors

        # Build transfer matrix in upright frame.
        cov_xx, cov_yy, cov_zz = cov_eig_vals
        RDx = scipy.special.elliprd(cov_yy, cov_zz, cov_xx)
        RDy = scipy.special.elliprd(cov_xx, cov_zz, cov_yy)
        RDz = scipy.special.elliprd(cov_xx, cov_yy, cov_zz)

        factor = 0.5 * self.sc_factor * ((1.0 / 5.0) ** 1.5)
        kappa_x = factor * RDx
        kappa_y = factor * RDy
        kappa_z = factor * RDz

        M = np.identity(7)
        M[1, 0] = kappa_x * length
        M[3, 2] = kappa_y * length
        M[5, 4] = kappa_z * length

        # Build matrix to undo x-y-z diagonalization.
        A = build_diag_matrix_from_xyz_eig(cov_eig_vecs)
        A_inv = A.T

        # Build matrix for translation to centroid.
        T = np.identity(7)
        for i in (0, 2, 4):
            T[i, -1] = centroid[i]

        T_inv = np.copy(T)
        T_inv[:-1, -1] = -T[:-1, -1]

        # Compute transfer matrix in lab frame.
        M = L @ T @ A @ M @ A_inv @ T_inv @ L_inv

        # Convert from z' to dE
        return convert_matrix_zp_to_dE(M, self.sync_part)


class EnvelopeTracker:
    def __init__(self, lattice: AccLattice, space_charge: str | None = None) -> None:
        self.lattice = lattice
        self.space_charge = space_charge

        for node in self.lattice.getNodes():
            if type(node) in (BendTEAPOT, BendLINAC):
                if node.getParam("ea1") != 0.0 or node.getParam("ea2") != 0.0:
                    message = f"Found bend ea1 or ea2 != 0.0 ({node.getName()}.)"
                    message += " Nonzero edge angles are not yet supported in envelope tracking."
                    message += " Setting ea1 and ea2 to 0.0."
                    warnings.warn(message)

                    node.setParam("ea1", 0.0)
                    node.setParam("ea2", 0.0)

    def track(self, envelope: Envelope) -> None:
        sync_part = envelope.sync_part
        charge = envelope.charge()

        for node_index, node in enumerate(self.lattice.getNodes()):
            for child_node in node.getChildNodes(ENTRANCE):
                matrix = track_sync_part(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    envelope.transform(matrix)

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=BEFORE):
                    matrix = track_sync_part(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        envelope.transform(matrix)

                matrix_sc = None
                if self.space_charge:
                    length = node.getLength(part_index)
                    if length > 0:
                        if self.space_charge == "2d":
                            matrix_sc = envelope.sc_matrix_2d(length)
                        elif self.space_charge == "3d":
                            matrix_sc = envelope.sc_matrix_3d(length)
                        else:
                            raise ValueError

                matrix = track_sync_part(node, sync_part=sync_part, charge=charge, index=part_index)
                if matrix is not None:
                    if matrix_sc is not None:
                        matrix = matrix @ matrix_sc
                    envelope.transform(matrix)

                for child_node in node.getChildNodes(BODY, part_index, place_in_part=AFTER):
                    matrix = track_sync_part(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        envelope.transform(matrix)

            for child_node in node.getChildNodes(EXIT):
                matrix = track_sync_part(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    envelope.transform(matrix)


    def track_history(self, envelope: Envelope) -> dict[str, list]:
        """Same as track but returns parameters vs. position in lattice."""
        history = {}
        history["position"] = []
        history["rms_x"] = []
        history["rms_y"] = []
        history["rms_z"] = []
        history["kin_energy"] = []

        sync_part = envelope.sync_part
        charge = envelope.charge()
        node_positions = self.lattice.getNodePositionsDict()

        history["position"].append(0.0)
        history["rms_x"].append(1000.0 * envelope.rms(0))
        history["rms_y"].append(1000.0 * envelope.rms(2))
        history["rms_z"].append(1000.0 * envelope.rms(4))
        history["kin_energy"].append(envelope.sync_part.kinEnergy())

        for node_index, node in enumerate(self.lattice.getNodes()):
            for child_node in node.getChildNodes(ENTRANCE):
                matrix = track_sync_part(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    envelope.transform(matrix)

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=BEFORE):
                    matrix = track_sync_part(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        envelope.transform(matrix)

                matrix_sc = None
                if self.space_charge:
                    length = node.getLength(part_index)
                    if length > 0:
                        if self.space_charge == "2d":
                            matrix_sc = envelope.sc_matrix_2d(length)
                        elif self.space_charge == "3d":
                            matrix_sc = envelope.sc_matrix_3d(length)
                        else:
                            raise ValueError

                matrix = track_sync_part(node, sync_part=sync_part, charge=charge, index=part_index)
                if matrix is not None:
                    if matrix_sc is not None:
                        matrix = matrix @ matrix_sc
                    envelope.transform(matrix)

                position_start, position_stop = node_positions[node]
                position = position_start + node.getLength(part_index) * (part_index + 1)

                history["position"].append(position)
                history["rms_x"].append(1000.0 * envelope.rms(0))
                history["rms_y"].append(1000.0 * envelope.rms(2))
                history["rms_z"].append(1000.0 * envelope.rms(4))
                history["kin_energy"].append(envelope.sync_part.kinEnergy())

                for child_node in node.getChildNodes(BODY, part_index, place_in_part=AFTER):
                    matrix = track_sync_part(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        envelope.transform(matrix)

            for child_node in node.getChildNodes(EXIT):
                matrix = track_sync_part(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    envelope.transform(matrix)

        return history