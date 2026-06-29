import numpy as np
import scipy.constants
import scipy.special

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle

from .matrix import convert_matrix_zp_to_dE
from .utils import gen_dist
from .utils import proj_cov_matrix


def build_diag_matrix_from_xyz_eig(eigenvectors: np.ndarray) -> np.ndarray:
    A = np.eye(7)
    for i in range(eigenvectors.shape[0]):
        for j in range(eigenvectors.shape[1]):
            row = i * 2
            col = j * 2
            A[row, col] = A[row + 1, col + 1] = eigenvectors[i, j]
    return A


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

        # [TO DO]
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

    def copy(self):
        return Envelope(
            bunch=self.bunch,
            cov_matrix=self.cov_matrix,
            centroid=self.centroid,
            intensity=self.intensity
        )

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

        # Get kick length in rest frame
        length_rest = length * gamma

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
        M[1, 0] = kappa_x * length_rest
        M[3, 2] = kappa_y * length_rest
        M[5, 4] = kappa_z * length_rest

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
