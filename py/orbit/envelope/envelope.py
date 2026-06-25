import math

import numpy as np
import scipy.special

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle
from orbit.lattice import AccNode
from orbit.lattice import AccLattice
from orbit.matrix_lattice.analytic import convert_matrix_zp_to_dE

from .utils import gen_dist
from .utils import proj_cov_matrix


ENTRANCE = AccNode.ENTRANCE
BODY = AccNode.BODY
EXIT = AccNode.EXIT

BEFORE = AccNode.BEFORE
AFTER = AccNode.AFTER

CLASSICAL_PROTON_RADIUS = 1.534697049469832e-18  # [m]


def build_diag_matrix_from_xyz_eig(eigenvectors: np.ndarray) -> np.ndarray:
    A = np.eye(7)
    for i in range(eigenvectors.shape[0]):
        for j in range(eigenvectors.shape[1]):
            row = i * 2
            col = j * 2
            A[row, col] = A[row + 1, col + 1] = eigenvectors[i, j]
    return A


class Envelope:
    """Represents beam envelope and centroid.

    Attributes:
        moment_matrix: 7 x 7 matrix containing first and second moments.
            Define the phase space vector X = [x, x', y, y', z, dE]^T and
            augmented vector Y = [x, x', y, y', z, dE, 1].

            Let X evolve according to X -> MX + U, where M is a 6 x 6 transfer matrix
            and U is 6 x 1 "driving" vector. The augmented vector Y evolves according
            to Y -> NY, where N = [[M, U], [0, 1]] is a 7 x 7 matrix.

            Let S = <YY^T> = [[<XX^T>, <X>], [<X^T>, 1]] = [[R, C], [C^T, 1]]. Here
            R = <XX^T> is the matrix of second moments, or "autocorrelation" matrix,
            and C = <X> is the mean/centroid vector. (To get the covariance matrix:
            <(X - C)(X - C)^T> = <XX^T> - <X><X>^T = R - CC^T.) S evolves according
            to S -> N S N^T.
        bunch: Bunch containing synchronous particle and (optionally) test particles.
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

        self.centroid = centroid
        if self.centroid is None:
            self.centroid = np.zeros(6)

        self.cov_matrix = cov_matrix
        if self.cov_matrix is None:
            self.cov_matrix = np.eye(6)

        self.intensity = 0.0
        self.set_intensity(intensity)

    def set_intensity(self, intensity: float) -> None:
        self.intensity = intensity
        self.sc_factor = (
            2.0
            * intensity
            * CLASSICAL_PROTON_RADIUS
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

    def transform(self, matrix: np.ndarray | None) -> None:
        if matrix is not None:
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

        bunch_length = 4.0 * np.sqrt(cov_matrix[4, 4])
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
        # z -> gamma * z
        # x' = dx/ds -> x' / gamma
        # y' = dy/ds -> y' / gamma
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
        kappa_x = factor * RDx  # [1 / m]
        kappa_y = factor * RDy  # [1 / m]
        kappa_z = factor * RDz  # [1 / m]

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

    def track(self, envelope: Envelope) -> None:
        charge = envelope.charge()
        for node_index, node in enumerate(self.lattice.getNodes()):
            for child_node in node.getChildNodes(ENTRANCE):
                matrix = child_node.matrix(sync_part=envelope.sync_part, charge=charge)
                envelope.transform(matrix)

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=BEFORE):
                    matrix = child_node.matrix(sync_part=envelope.sync_part, charge=charge)
                    envelope.transform(matrix)

                if self.space_charge:
                    length = node.getLength(part_index)
                    if self.space_charge == "2d":
                        matrix = envelope.sc_matrix_2d(length)
                    elif self.space_charge == "3d":
                        matrix = envelope.sc_matrix_3d(length)
                    else:
                        raise ValueError(f"Invalid space charge model: {self.space_charge}")
                    envelope.transform(matrix)

                matrix = node.matrix(sync_part=envelope.sync_part, charge=charge, index=part_index)
                envelope.transform(matrix)

                for child_node in node.getChildNodes(BODY, part_index, place_in_part=AFTER):
                    matrix = child_node.matrix(sync_part=envelope.sync_part, charge=charge)
                    envelope.transform(matrix)

            for child_node in node.getChildNodes(EXIT):
                matrix = child_node.matrix(sync_part=envelope.sync_part, charge=charge)
                envelope.transform(matrix)

    def track_history(self, envelope: Envelope) -> dict[str, list]:
        history = {}
        history["position"] = []
        history["rms_x"] = []
        history["rms_y"] = []
        history["rms_z"] = []

        charge = envelope.charge()
        node_positions = self.lattice.getNodePositionsDict()

        history["position"].append(0.0)
        history["rms_x"].append(1000.0 * envelope.rms(0))
        history["rms_y"].append(1000.0 * envelope.rms(2))
        history["rms_z"].append(1000.0 * envelope.rms(4))
        
        for node_index, node in enumerate(self.lattice.getNodes()):
            for child_node in node.getChildNodes(ENTRANCE):
                matrix = child_node.matrix(sync_part=envelope.sync_part, charge=charge)
                envelope.transform(matrix)

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=BEFORE):
                    matrix = child_node.matrix(sync_part=envelope.sync_part, charge=charge)
                    envelope.transform(matrix)

                if self.space_charge:
                    length = node.getLength(part_index)
                    if self.space_charge == "2d":
                        matrix = envelope.sc_matrix_2d(length)
                    elif self.space_charge == "3d":
                        matrix = envelope.sc_matrix_3d(length)
                    else:
                        raise ValueError(f"Invalid space charge model: {self.space_charge}")
                    envelope.transform(matrix)

                matrix = node.matrix(sync_part=envelope.sync_part, charge=charge, index=part_index)
                envelope.transform(matrix)

                position_start, position_stop = node_positions[node]
                position = position_start + node.getLength(part_index)

                history["position"].append(position)
                history["rms_x"].append(1000.0 * envelope.rms(0))
                history["rms_y"].append(1000.0 * envelope.rms(2))
                history["rms_z"].append(1000.0 * envelope.rms(4))

                for child_node in node.getChildNodes(BODY, part_index, place_in_part=AFTER):
                    matrix = child_node.matrix(sync_part=envelope.sync_part, charge=charge)
                    envelope.transform(matrix)

            for child_node in node.getChildNodes(EXIT):
                matrix = child_node.matrix(sync_part=envelope.sync_part, charge=charge)
                envelope.transform(matrix)

        return history