import math

import numpy as np

from ..core.bunch import SyncParticle
from ..lattice import AccNode
from ..lattice import AccLattice

from .matrix import MatrixFactory
from .utils import gen_dist


ENTRANCE = AccNode.ENTRANCE
BODY = AccNode.BODY
EXIT = AccNode.EXIT

BEFORE = AccNode.BEFORE
AFTER = AccNode.AFTER


def get_perveance(mass: float, kin_energy: float, line_density: float) -> float:
    classical_proton_radius = 1.53469e-18  # [m]
    gamma = 1.0 + (kin_energy / mass)  # Lorentz factor
    beta = np.sqrt(1.0 - (1.0 / gamma) ** 2)  # velocity/speed_of_light
    return (classical_proton_radius * line_density) / (beta**2 * gamma**3)


class Envelope:
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
        self.perveance = 0.0
        self.set_intensity(intensity)

    def set_intensity(self, intensity: float) -> None:
        self.intensity = intensity
        cov_matrix = self.cov()
        length = 2.0 * math.sqrt(cov_matrix[4, 4])  # assume uniform density
        self.perveance = get_perveance(
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


class EnvelopeTracker:
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
                    cov_matrix = envelope.cov()

                    if self.space_charge == "2d":
                        matrix = self.matrix_factory.space_charge_2d(
                            length=length, cov_matrix=cov_matrix, perveance=envelope.perveance
                        )
                    elif self.space_charge == "3d":
                        matrix = self.matrix_factory.space_charge_3d(
                            length=length, cov_matrix=cov_matrix, intensity=envelope.intensity
                        )
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
