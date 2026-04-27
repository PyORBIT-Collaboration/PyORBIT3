import math

import numpy as np

from ..core.bunch import SyncParticle
from ..lattice import AccNode
from ..lattice import AccLattice

from .matrix import MatrixFactory
from .utils import gen_dist
from .utils import proj_cov_matrix


ENTRANCE = AccNode.ENTRANCE
BODY = AccNode.BODY
EXIT = AccNode.EXIT

BEFORE = AccNode.BEFORE
AFTER = AccNode.AFTER


def get_perveance_2d(mass: float, kin_energy: float, line_density: float) -> float:
    classical_proton_radius = 1.534697049469832e-18  # [m]
    gamma = 1.0 + (kin_energy / mass)  # Lorentz factor
    beta = np.sqrt(1.0 - (1.0 / gamma) ** 2)  # velocity/speed_of_light
    return (2.0 * classical_proton_radius * line_density) / (beta**2 * gamma**3)


def get_perveance_3d():
    raise NotImplementedError()


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
        self.perveance_3d = None
        self.set_intensity(intensity)

    def set_intensity(self, intensity: float) -> None:
        self.intensity = intensity
        cov_matrix = self.cov()
        length = 4.0 * math.sqrt(cov_matrix[4, 4])  # assume uniform density
        self.perveance_2d = get_perveance_2d(
            mass=self.sync_part.mass(),
            kin_energy=self.sync_part.kinEnergy(),
            line_density=(self.intensity / length),
        )

    def centroid(self) -> np.ndarray:
        return self.matrix[0:6, 6]

    def cov(self) -> np.ndarray:
        return self.matrix[0:6, 0:6]

    def rms(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov()))

    def apply_transfer_matrix(self, transfer_matrix: np.ndarray) -> None:
        self.matrix = np.linalg.multi_dot([transfer_matrix, self.matrix, transfer_matrix.T])

    def sample(self, n: int, dist: str = "kv") -> np.ndarray:
        # Issue: covariance matrix is becoming non semi-positive definite,
        # giving error in cholesky decomposition.
        particles = gen_dist(n=n, cov_matrix=self.cov(), name=dist)
        particles = particles + self.centroid()
        return particles
    
    def sc_transfer_matrix_2d(self, length: float) -> np.ndarray:
        centroid = self.centroid()

        cov_matrix = self.cov()
        cov_matrix_proj = proj_cov_matrix(cov_matrix, axis=(0, 2))
        
        eig_res = np.linalg.eig(cov_matrix_proj)

        rx, ry = 2.0 * np.sqrt(eig_res.eigenvalues)      

        kappa_x = 2.0 * self.perveance_2d / (rx * (rx + ry))
        kappa_y = 2.0 * self.perveance_2d / (ry * (rx + ry))
        
        M = np.identity(7)
        M[1, 0] = kappa_x * length
        M[3, 2] = kappa_y * length

        A = build_diag_matrix_from_xyz_eig(eig_res.eigenvectors)
        
        T = np.identity(7)
        T[0, -1] = centroid[0]
        T[2, -1] = centroid[2] 

        V = np.matmul(T, A)
        V_inv = np.linalg.inv(V)
        return np.linalg.multi_dot([V, M, V_inv])
    
    def sc_transfer_matrix_3d(self, length: float) -> np.ndarray:
        centroid = self.centroid()

        cov_matrix = self.cov()
        cov_matrix_proj = proj_cov_matrix(cov_matrix, axis=(0, 2, 4))
        
        eig_res = np.linalg.eig(cov_matrix_proj)

        cov_xx, cov_yy, cov_zz = np.sqrt(eig_res.eigenvalues)      

        kappa_x = ... 
        kappa_y = ... 
        kappa_z = ...
        
        M = np.identity(7)
        M[1, 0] = kappa_x * length
        M[3, 2] = kappa_y * length
        M[5, 4] = kappa_z * length

        A = build_diag_matrix_from_xyz_eig(eig_res.eigenvectors)
        
        T = np.identity(7)
        T[0, -1] = centroid[0]
        T[2, -1] = centroid[2] 
        T[4, -1] = centroid[4]

        V = np.matmul(T, A)
        V_inv = np.linalg.inv(V)
        return np.linalg.multi_dot([V, M, V_inv])


class EnvelopeTracker:
    """Tracks envelope through linear lattice with optional linear space charge kicks."""
    def __init__(self, lattice: AccLattice, space_charge: str | None = None) -> None:
        self.lattice = lattice
        self.matrix_factory = MatrixFactory()
        self.space_charge = space_charge

    def track(self, envelope: Envelope) -> None:
        for node in self.lattice.getNodes():
            # Child nodes before node
            for child_node in node.getChildNodes(ENTRANCE):
                envelope.apply_transfer_matrix(self.matrix_factory(child_node, envelope.sync_part))

            for part_index in range(node.getnParts()):
                # Child nodes before part
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=BEFORE):
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
                        raise ValueError(f"Invalid space charge model: {self.space_charge}")

                    envelope.apply_transfer_matrix(matrix)

                # Main node part
                envelope.apply_transfer_matrix(
                    self.matrix_factory(node, envelope.sync_part, part_index)
                )

                # Child nodes after part
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=AFTER):
                    envelope.apply_transfer_matrix(
                        self.matrix_factory(child_node, envelope.sync_part)
                    )

            # Child nodes after node
            for child_node in node.getChildNodes(EXIT):
                envelope.apply_transfer_matrix(self.matrix_factory(child_node, envelope.sync_part))
