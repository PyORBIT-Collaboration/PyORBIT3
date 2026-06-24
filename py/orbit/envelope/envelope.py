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
        moment_matrix: 7 x 7 covariance matrix for augmented phase space vector.
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

        self.sync_part = bunch.getSyncParticle()
        self.dim = 6

        if centroid is None:
            centroid = np.zeros(6)

        if cov_matrix is None:
            cov_matrix = np.eye(6)

        self.moment_matrix = np.zeros((7, 7))
        self.moment_matrix[: self.dim, : self.dim] = cov_matrix + np.outer(
            centroid, centroid
        )
        self.moment_matrix[: self.dim, self.dim] = centroid
        self.moment_matrix[self.dim, : self.dim] = centroid
        self.moment_matrix[self.dim, self.dim] = 1.0

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

    def centroid(self) -> np.ndarray:
        return np.copy(self.moment_matrix[: self.dim, self.dim])

    def cov(self) -> np.ndarray:
        autocorrelation_matrix = self.moment_matrix[: self.dim, : self.dim]
        centroid = self.moment_matrix[: self.dim, self.dim]
        return autocorrelation_matrix - np.outer(centroid, centroid)

    def rms(self, axis: int = None) -> float | np.ndarray:
        rms_arr = np.sqrt(np.diag(self.cov()))
        return rms_arr[axis]

    def apply_transfer_matrix(self, transfer_matrix: np.ndarray | None) -> None:
        if transfer_matrix is not None:
            self.moment_matrix = (
                transfer_matrix @ self.moment_matrix @ transfer_matrix.T
            )

    def sample(self, size: int, dist: str = "kv") -> np.ndarray:
        # Issue: covariance matrix is becoming non semi-positive definite,
        # giving error in cholesky decomposition.
        particles = gen_dist(size=size, cov_matrix=self.cov(), name=dist)
        particles = particles + self.centroid()
        return particles

    def sc_transfer_matrix_2d(self, length: float) -> np.ndarray:
        # Extract beam centroid and covariance matrix.
        centroid = self.centroid()
        cov_matrix = self.cov()

        # Project covariance matrix onto x-y plane.
        cov_matrix_proj = proj_cov_matrix(cov_matrix, axis=(0, 2))

        # Compute eigenvalues and eigenvectors of x-y covariance matrix.
        cov_eig_res = np.linalg.eigh(cov_matrix_proj)
        cov_eig_vals = cov_eig_res.eigenvalues
        cov_eig_vecs = cov_eig_res.eigenvectors

        # Compute rms beam sizes in upright frame.
        rx = 2.0 * math.sqrt(cov_eig_vals[0])
        ry = 2.0 * math.sqrt(cov_eig_vals[1])

        # Build transfer matrix in upright frame.
        bunch_length = 4.0 * self.rms(axis=4)
        perveance = self.sc_factor / bunch_length
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

        # Compute transfer matrix in lab frame.
        V = np.matmul(T, A)
        V_inv = np.linalg.inv(V)
        return np.linalg.multi_dot([V, M, V_inv])

    def sc_transfer_matrix_3d(self, length: float) -> np.ndarray:
        # Build Lorentz matrix: rest frame to lab frame.
        # x -> x
        # y -> y
        # z -> gamma * z
        # x' = dx/ds -> x' / gamma
        # y' = dy/ds -> y' / gamma
        # z' = dz/ds -> z'
        lorentz_matrix = np.identity(7)
        lorentz_matrix[1, 1] = self.gamma()
        lorentz_matrix[3, 3] = self.gamma()
        lorentz_matrix[4, 4] = 1.0 / self.gamma()
        lorentz_matrix_inv = np.linalg.inv(lorentz_matrix)

        # Get centroid in rest frame.
        centroid = self.centroid()
        centroid[4] *= self.gamma()

        # Get covariance matrix in rest frame.
        cov_matrix = self.cov()
        cov_matrix = np.linalg.multi_dot(
            [lorentz_matrix_inv[:-1, :-1], cov_matrix, lorentz_matrix_inv[:-1, :-1].T]
        )

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

        # Build matrix for translation to centroid.
        T = np.identity(7)
        for i in (0, 2, 4):
            T[i, -1] = centroid[i]

        # Build matrix for Lorentz boost (length contraction).
        L = lorentz_matrix

        # Compute transfer matrix in lab frame.
        V = np.linalg.multi_dot([L, T, A])
        M = np.linalg.multi_dot([V, M, np.linalg.inv(V)])

        # Convert from z' to dE
        return convert_matrix_zp_to_dE(M, self.sync_part)


class EnvelopeTracker:
    def __init__(self, lattice: AccLattice, space_charge: str | None = None) -> None:
        self.lattice = lattice
        self.space_charge = space_charge

    def track(self, envelope: Envelope) -> None:
        for node in self.lattice.getNodes():
            for child_node in node.getChildNodes(ENTRANCE):
                envelope.apply_transfer_matrix(child_node.matrix(envelope.sync_part))

            for index in range(node.getnParts()):
                for child_node in node.getChildNodes(BODY, index, place_in_part=BEFORE):
                    envelope.apply_transfer_matrix(
                        child_node.matrix(envelope.sync_part)
                    )

                if self.space_charge:
                    length = node.getLength(index)
                    if self.space_charge == "2d":
                        matrix = envelope.sc_transfer_matrix_2d(length)
                    elif self.space_charge == "3d":
                        matrix = envelope.sc_transfer_matrix_3d(length)
                    else:
                        raise ValueError(
                            f"Invalid space charge model: {self.space_charge}"
                        )
                    envelope.apply_transfer_matrix(matrix)

                envelope.apply_transfer_matrix(node.matrix(envelope.sync_part, index))

                for child_node in node.getChildNodes(BODY, index, place_in_part=AFTER):
                    envelope.apply_transfer_matrix(
                        child_node.matrix(envelope.sync_part)
                    )

            for child_node in node.getChildNodes(EXIT):
                envelope.apply_transfer_matrix(child_node.matrix(envelope.sync_part))
