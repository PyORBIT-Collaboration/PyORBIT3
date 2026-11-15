import copy
import math

import numpy as np
import scipy.constants
import scipy.special

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis

from .matrix import get_zp_coeff
from .utils import get_classical_radius
from .utils import sample_distribution


# For a uniform one-dimensional distribution over length L, the standard
# deviation is L / sqrt(12). This quantity is used to calculate the line
# density for two-dimensional space charge kicks.
RMS_BUNCH_LENGTH_FACTOR = math.sqrt(12.0)


def get_bunch_covariance(bunch: Bunch) -> np.ndarray:
    twiss_calc = BunchTwissAnalysis()
    twiss_calc.analyzeBunch(bunch)
    cov_matrix = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            cov_matrix[i, j] = cov_matrix[j, i] = twiss_calc.getCorrelation(i, j)
    return cov_matrix


def get_bunch_centroid(bunch: Bunch) -> np.ndarray:
    twiss_calc = BunchTwissAnalysis()
    twiss_calc.analyzeBunch(bunch)
    centroid = np.zeros(6)
    for i in range(6):
        centroid[i] = twiss_calc.getAverage(i)
    return centroid


class Envelope:
    """Represents the first and second moments of a phase-space distribution.."""

    def __init__(
        self,
        cov_matrix: np.ndarray | None = None,
        centroid: np.ndarray | None = None,
        kin_energy: float | None = None,
        mass: float | None = None,
        charge: float = 1.0,
        time: float = 0.0,
        intensity: float = 0.0,
        bunch: Bunch | None = None,
    ) -> None:
        """Constructor.

        Args:
            cov_matrix: 6 x 6 covariance matrix.
            centroid: 6 x 1 centroid vector.
            kin_energy: Synchronous particle kinetic energy [GeV].
            mass: Particle mass [GeV / c^2].
            charge: Particle charge [e].
            time: Synchronous particle time [s].
            intensity: Total number of particles in bunch.
            bunch: If provided, the parameters above are calculated
                from the bunch parameters and macroparticle coordiantes.
        """
        if bunch is not None:
            cov_matrix = get_bunch_covariance(bunch)
            centroid = get_bunch_centroid(bunch)
            mass = bunch.mass()
            charge = bunch.charge()
            sync_part = bunch.getSyncParticle()
            kin_energy = sync_part.kinEnergy()
            time = sync_part.time()
            intensity = bunch.macroSize() * bunch.getSizeGlobal()

        self.centroid = np.zeros(6) if centroid is None else centroid
        self.cov_matrix = np.eye(6) if cov_matrix is None else cov_matrix
        self.kin_energy = kin_energy
        self.mass = mass
        self.charge = charge
        self.time = time
        self.intensity = intensity
        self.classical_radius = get_classical_radius(self.charge, self.mass)
        self.rms_bunch_length_factor = RMS_BUNCH_LENGTH_FACTOR

    def copy(self) -> "Envelope":
        return copy.deepcopy(self)

    @property
    def gamma(self) -> float:
        return 1.0 + self.kin_energy / self.mass

    @property
    def beta(self) -> float:
        return self.momentum / (self.kin_energy + self.mass)

    @property
    def momentum(self) -> float:
        return math.sqrt(self.kin_energy * (self.kin_energy + 2.0 * self.mass))

    @property
    def rigidity(self) -> float:
        return 3.335640952 * self.momentum / self.charge

    @property
    def sc_factor(self) -> float:
        return 2.0 * self.intensity * self.classical_radius / (self.beta**2 * self.gamma**3)

    def rms(self, axis: int | None = None) -> float | np.ndarray:
        return np.sqrt(np.diag(self.cov_matrix))[axis]

    def apply_matrix(self, matrix: np.ndarray) -> None:
        m = matrix[:-1, :-1]
        u = matrix[:-1, -1]
        self.cov_matrix = m @ self.cov_matrix @ m.T
        self.centroid = m @ self.centroid + u

    def sample(self, size: int, dist: str = "kv") -> np.ndarray:
        particles = sample_distribution(size=size, cov_matrix=self.cov_matrix, name=dist)
        particles = particles + self.centroid
        return particles

    def to_bunch(self, size: int, dist: str = "gauss") -> Bunch:
        bunch = Bunch()
        bunch.mass(self.mass)
        bunch.charge(self.charge)
        bunch.getSyncParticle().kinEnergy(self.kin_energy)
        bunch.getSyncParticle().time(self.time)
        bunch.macroSize(self.intensity / size)
        if size > 0:
            particles = self.sample(size=size, dist=dist)
            for i in range(particles.shape[0]):
                bunch.addParticle(*particles[i])
        return bunch

    def sc_matrix_2d(self, length: float) -> np.ndarray:
        cov_matrix = self.cov_matrix
        cov_xx = cov_matrix[0, 0]
        cov_yy = cov_matrix[2, 2]
        cov_xy = cov_matrix[0, 2]

        # Calculate kick strengths in upright frame.
        phi = -0.5 * math.atan2(2.0 * cov_xy, cov_xx - cov_yy)
        sin_phi = math.sin(phi)
        cos_phi = math.cos(phi)
        rx = 2.0 * math.sqrt(
            abs(cov_xx * cos_phi**2 + cov_yy * sin_phi**2 - 2.0 * cov_xy * sin_phi * cos_phi)
        )
        ry = 2.0 * math.sqrt(
            abs(cov_xx * sin_phi**2 + cov_yy * cos_phi**2 + 2.0 * cov_xy * sin_phi * cos_phi)
        )

        bunch_length = self.rms_bunch_length_factor * math.sqrt(cov_matrix[4, 4])
        factor = 2.0 * self.sc_factor * length / (bunch_length * (rx + ry))
        kx = factor / rx
        ky = factor / ry

        # Rotate to lab frame.
        kxx = kx * cos_phi**2 + ky * sin_phi**2
        kxy = (ky - kx) * sin_phi * cos_phi
        kyy = kx * sin_phi**2 + ky * cos_phi**2
        matrix = np.eye(7)
        matrix[1, 0] = kxx
        matrix[1, 2] = matrix[3, 0] = kxy
        matrix[3, 2] = kyy

        # Account for nonzero centroid.
        x = self.centroid[0]
        y = self.centroid[2]
        matrix[1, 6] = -kxx * x - kxy * y
        matrix[3, 6] = -kxy * x - kyy * y
        return matrix

    def sc_matrix_3d(self, length: float) -> np.ndarray:
        gamma = self.gamma

        # Get position covariance with z boosted to rest frame.
        cov_position = self.cov_matrix[::2, ::2].copy()
        cov_position[2, :] *= gamma
        cov_position[:, 2] *= gamma
        eigenvalues, eigenvectors = np.linalg.eigh(cov_position)
        cov_xx, cov_yy, cov_zz = eigenvalues

        strengths = np.array(
            [
                scipy.special.elliprd(cov_yy, cov_zz, cov_xx),
                scipy.special.elliprd(cov_xx, cov_zz, cov_yy),
                scipy.special.elliprd(cov_xx, cov_yy, cov_zz),
            ]
        )
        strengths *= 0.5 * self.sc_factor * (0.2**1.5) * length

        # Rotate out of upright frame.
        kick = (eigenvectors * strengths) @ eigenvectors.T

        # Transform to lab frame and convert z' to dE.
        kick[:2, :] *= gamma
        kick[:, 2] *= gamma
        kick[2, :] /= get_zp_coeff(self)

        matrix = np.eye(7)
        matrix[1:6:2, :6:2] = kick
        matrix[1:6:2, 6] = -kick @ self.centroid[::2]
        return matrix

    def sc_matrix(self, length: float, model: str) -> np.ndarray:
        if model == "2d":
            return self.sc_matrix_2d(length)
        elif model == "3d":
            return self.sc_matrix_3d(length)
        else:
            raise ValueError(f"Unknown space charge model {model}")
