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


def calc_perveance_2d(gamma: float, line_density: float) -> float:
    classical_proton_radius = 1.534697049469832e-18  # [m]
    beta = np.sqrt(1.0 - (1.0 / gamma) ** 2)
    return (2.0 * classical_proton_radius * line_density) / (beta**2 * gamma**3)


def calc_perveance_3d(gamma: float, mass: float, total_charge: float) -> float:
    bg2 = gamma**2 - 1.0
    return total_charge * (1.0e-7 * speed_of_light**2) * (1.0 / (gamma * bg2)) * abs(charge_electron) / mass


def build_diag_matrix_from_xyz_eig(eigenvectors: np.ndarray) -> np.ndarray:
    A = np.eye(7)
    for i in range(eigenvectors.shape[0]):
        for j in range(eigenvectors.shape[1]):
            row = i * 2
            col = j * 2
            A[row, col] = A[row + 1, col + 1] = eigenvectors[i, j]
    return A


def build_sc_matrix_2d(
    cov_xx: float,
    cov_yy: float,
    perveance: float,
    length: float,
) -> np.ndarray:
    """Return 7 x 7 space charge matrix for upright 2D ellipsoid."""
    r_x = 2.0 * math.sqrt(cov_xx)
    r_y = 2.0 * math.sqrt(cov_yy)
    kappa_x = 2.0 * perveance / (r_x * (r_x + r_y))
    kappa_y = 2.0 * perveance / (r_y * (r_x + r_y))
    
    matrix = np.identity(7)
    matrix[1, 0] = kappa_x * length
    matrix[3, 2] = kappa_y * length
    return matrix


def build_sc_matrix_3d(
    cov_xx: float, 
    cov_yy: float, 
    cov_zz: float, 
    perveance: float, 
    length: float,
    gamma: float, 
) -> np.ndarray:    
    """Return 7 x 7 space charge matrix for upright 3D ellipsoid."""

    cov_xx = cov_xx
    cov_yy = cov_yy
    cov_zz = cov_zz * gamma * gamma
    
    scale_factor = 5.0 ** 1.5
    RDx = scipy.special.elliprd(cov_yy, cov_zz, cov_xx) / scale_factor
    RDy = scipy.special.elliprd(cov_zz, cov_xx, cov_yy) / scale_factor
    RDz = scipy.special.elliprd(cov_xx, cov_yy, cov_zz) / scale_factor
  
    kx = gamma * length * perveance * RDx
    ky = gamma * length * perveance * RDy
    kz = gamma * length * perveance * RDz

    matrix = np.identity(7)
    matrix[1, 0] = kx
    matrix[3, 2] = ky
    matrix[5, 4] = kz
    return matrix


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
        cov_matrix = self.cov()
        length = 4.0 * math.sqrt(cov_matrix[4, 4])  # assume uniform density
        self.perveance_2d = calc_perveance_2d(
            gamma=self.gamma(),
            line_density=(self.intensity / length),
        )
        self.perveance_3d = calc_perveance_3d(
            gamma=self.gamma(),
            mass=self.mass(),
            total_charge=(self.intensity * charge_electron),
        )

    def gamma(self) -> float:
        return self.sync_part.gamma()
    
    def mass(self) -> float:
        return self.sync_part.mass()

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

        M = build_sc_matrix_2d(
            cov_xx=eig_res.eigenvalues[0], 
            cov_yy=eig_res.eigenvalues[1], 
            perveance=self.perveance_2d, 
            length=length
        )

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

        M = build_sc_matrix_3d(
            cov_xx=eig_res.eigenvalues[0], 
            cov_yy=eig_res.eigenvalues[1], 
            cov_zz=eig_res.eigenvalues[2],
            perveance=self.perveance_3d, 
            length=length,
            gamma=self.gamma(),
        )   

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
