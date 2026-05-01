import math

import numpy as np
import scipy.special

from ..core.bunch import SyncParticle
from ..lattice import AccNode
from ..lattice import AccLattice
from ..utils.consts import speed_of_light
from ..utils.consts import charge_electron

from .matrix import MatrixFactory
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
        matrix: 7 x 7 covariance matrix for augmented phase space vector.

            Define the phase space vector X = [x, x', y, y', z, dE]^T and
            augmented vector Y = [x, x', y, y', z, dE, 1].

            Let X evolve according to X -> MX + U, where M is a 6 x 6 transfer matrix
            and U is 6 x 1 "driving" vector. The augmented vector Y evolves according
            to Y -> NY, where N = [[M, U], [0, 1]] is a 7 x 7 matrix.

            We track the 7 x 7 covariance matrix of Y:

            R = <YY^T> = [[<XX^T>, <X>], [<X^T>, 1]],

            which contains both the phase space covariance matrix and centroid vector.
            R evolves according to R -> N R N^T.
    """

    def __init__(
        self,
        sync_part: SyncParticle,
        cov_matrix: np.ndarray = None,
        centroid: np.ndarray = None,
        intensity: float = 0.0,
    ) -> None:
        self.sync_part = sync_part

        if centroid is None:
            centroid = np.zeros(6)

        if cov_matrix is None:
            cov_matrix = np.eye(6)

        self.matrix = np.zeros((7, 7))
        self.matrix[0:6, 0:6] = cov_matrix
        self.matrix[0:6, 6] = centroid
        self.matrix[6, 0:6] = centroid
        self.matrix[6, 6] = 1.0

        self.intensity = 0.0
        self.perveance_2d = 0.0
        self.perveance_3d = 0.0
        self.set_intensity(intensity)

    def set_intensity(self, intensity: float) -> None:
        self.intensity = intensity
        self.perveance = (
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

    def centroid(self) -> np.ndarray:
        return np.copy(self.matrix[0:6, 6])

    def cov(self) -> np.ndarray:
        return np.copy(self.matrix[0:6, 0:6])

    def rms(self, axis: int = None) -> float | np.ndarray:
        rms_arr = np.sqrt(np.diag(self.cov()))
        return rms_arr[axis]

    def apply_transfer_matrix(self, transfer_matrix: np.ndarray) -> None:
        self.matrix = np.linalg.multi_dot(
            [transfer_matrix, self.matrix, transfer_matrix.T]
        )

    def sample(self, n: int, dist: str = "kv") -> np.ndarray:
        # Issue: covariance matrix is becoming non semi-positive definite,
        # giving error in cholesky decomposition.
        particles = gen_dist(n=n, cov_matrix=self.cov(), name=dist)
        particles = particles + self.centroid()
        return particles

    def sc_transfer_matrix_2d(self, length: float) -> np.ndarray:
        # Extract beam centroid and covariance matrix.
        centroid = self.centroid()
        cov_matrix = self.cov()

        # Project covariance matrix onto x-y plane.
        cov_matrix_proj = proj_cov_matrix(cov_matrix, axis=(0, 2))

        # Compute eigenvalues and eigenvectors of x-y covariance matrix.
        cov_eig_res = np.linalg.eig(cov_matrix_proj)
        cov_eig_vals = cov_eig_res.eigenvalues
        cov_eig_vecs = cov_eig_res.eigenvectors

        # Compute rms beam sizes in upright frame.
        rx = 2.0 * math.sqrt(cov_eig_vals[0])
        ry = 2.0 * math.sqrt(cov_eig_vals[1])

        # Build transfer matrix in upright frame.
        bunch_length = 4.0 * self.rms(axis=4)
        perveance = self.perveance / bunch_length
        factor = 2.0 * perveance / (rx + ry)
        kappa_x = factor / rx
        kappa_y = factor / ry

        M = np.identity(7)
        M[1, 0] = kappa_x * length
        M[3, 2] = kappa_y * length

        # Build matrix to undo x-y diagonalization.
        A = build_diag_matrix_from_xyz_eig(cov_eig_vecs)

        # Build matrix to translate to centroid.
        T = np.identity(7)
        T[0, -1] = centroid[0]
        T[2, -1] = centroid[2]

        # Compute matrix in lab frame.
        #   x = V u = T A u.
        #   u -> M u
        #   x -> V M V^-1 x
        V = np.matmul(T, A)
        V_inv = np.linalg.inv(V)
        return np.linalg.multi_dot([V, M, V_inv])

    def sc_transfer_matrix_3d(self, length: float) -> np.ndarray:
        centroid = self.centroid()
        centroid[4] *= self.gamma()

        cov_matrix = self.cov()
        cov_matrix[4, 4] *= self.gamma() ** 2
        cov_matrix_proj = proj_cov_matrix(cov_matrix, axis=(0, 2, 4))

        # Compute eigenvalues and eigenvectors of x-y covariance matrix.
        cov_eig_res = np.linalg.eig(cov_matrix_proj)
        cov_eig_vals = cov_eig_res.eigenvalues
        cov_eig_vecs = cov_eig_res.eigenvectors

        # Build transfer matrix in upright frame.
        cov_xx, cov_yy, cov_zz = cov_eig_vals
        RDx = scipy.special.elliprd(cov_yy, cov_zz, cov_xx)
        RDy = scipy.special.elliprd(cov_xx, cov_zz, cov_yy)
        RDz = scipy.special.elliprd(cov_xx, cov_yy, cov_zz)

        factor = 0.5 * self.perveance * ((1.0 / 5.0) ** 1.5)
        kappa_x = factor * RDx  # [1 / m]
        kappa_y = factor * RDy  # [1 / m]
        kappa_z = factor * RDz  # [1 / m]
        kappa_z *= self.gamma() ** 3 * self.beta() ** 2 * self.mass()  # [GeV / m]

        M = np.identity(7)
        M[1, 0] = kappa_x * length
        M[3, 2] = kappa_y * length
        M[5, 4] = kappa_z * length

        # Build matrix to undo x-y-z diagonalization.
        A = build_diag_matrix_from_xyz_eig(cov_eig_vecs)

        # Build matrix for translation to centroid.
        T = np.identity(7)
        for i in (0, 2, 4):
            T[i, -1] = centroid[i]

        # Build matrix for Lorentz boost (length contraction).
        L = np.identity(7)
        L[4, 4] = 1.0 / self.gamma()

        # Compute matrix in lab frame.
        #   x = L V u = L T A u.
        #   u -> M u
        #   x -> V M V^-1 x
        V = np.matmul(T, A)
        V_inv = np.linalg.inv(V)
        return np.linalg.multi_dot([V, M, V_inv])


class EnvelopeTracker:
    """Tracks envelope through linear lattice with optional linear space charge kicks."""

    def __init__(self, lattice: AccLattice, space_charge: str | None = None, ignore_unknown: bool = False) -> None:
        self.lattice = lattice
        self.matrix_factory = MatrixFactory(ignore_unknown=ignore_unknown)
        self.space_charge = space_charge

    def track(self, envelope: Envelope) -> None:
        for node in self.lattice.getNodes():
            # Child nodes before node
            for child_node in node.getChildNodes(ENTRANCE):
                envelope.apply_transfer_matrix(
                    self.matrix_factory(child_node, envelope.sync_part)
                )

            for part_index in range(node.getnParts()):
                # Child nodes before part
                for child_node in node.getChildNodes(
                    BODY, part_index, place_in_part=BEFORE
                ):
                    envelope.apply_transfer_matrix(
                        self.matrix_factory(child_node, envelope.sync_part)
                    )

                # Space charge
                if self.space_charge:
                    length = node.getLength(part_index)
                    if self.space_charge == "2d":
                        matrix = envelope.sc_transfer_matrix_2d(length)
                    elif self.space_charge == "3d":
                        matrix = envelope.sc_transfer_matrix_3d(length)
                    else:
                        raise ValueError(
                            f"Invalid space charge model: {self.space_charge}"
                        )

                    # print("debug space charge matrix")
                    # print(matrix)

                    envelope.apply_transfer_matrix(matrix)

                # Main node part
                envelope.apply_transfer_matrix(
                    self.matrix_factory(node, envelope.sync_part, part_index)
                )

                # Child nodes after part
                for child_node in node.getChildNodes(
                    BODY, part_index, place_in_part=AFTER
                ):
                    envelope.apply_transfer_matrix(
                        self.matrix_factory(child_node, envelope.sync_part)
                    )

            # Child nodes after node
            for child_node in node.getChildNodes(EXIT):
                envelope.apply_transfer_matrix(
                    self.matrix_factory(child_node, envelope.sync_part)
                )
