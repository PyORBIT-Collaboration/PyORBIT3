import math

import numpy as np

from ..core.bunch import SyncParticle
from ..lattice import AccNode
from ..lattice import AccLattice

from .matrix import MatrixFactory


ENTRANCE = AccNode.ENTRANCE
BODY = AccNode.BODY
EXIT = AccNode.EXIT

BEFORE = AccNode.BEFORE
AFTER = AccNode.AFTER


class Envelope:
    def __init__(self, sync_part: SyncParticle, cov_matrix: np.ndarray = None, centroid: np.ndarray = None) -> None:
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

    def centroid(self) -> np.ndarray:
        return self.matrix[0:6, 6]
    
    def cov(self) -> np.ndarray:
        return self.matrix[0:6, 0:6]
    
    def rms(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov()))

    def apply_transfer_matrix(self, transfer_matrix: np.ndarray) -> None:
        self.matrix = np.linalg.multi_dot([transfer_matrix, self.matrix, transfer_matrix.T])

    def space_charge_matrix(self, length: float) -> np.ndarray:
        """Return transfer matrix from linear space charge kick."""
        raise NotImplementedError()
    

class EnvelopeTracker:
    def __init__(self, lattice: AccLattice) -> None:
        self.lattice = lattice
        self.matrix_factory = MatrixFactory()

    def track(self, envelope: Envelope) -> None:
        for node in self.lattice.getNodes():
            for child_node in node.getChildNodes(ENTRANCE):
                envelope.apply_transfer_matrix(
                    self.matrix_factory(child_node, envelope.sync_part)
                )

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=BEFORE):
                    envelope.apply_transfer_matrix(
                        self.matrix_factory(child_node, envelope.sync_part)
                    )

                envelope.apply_transfer_matrix(
                    self.matrix_factory(node, envelope.sync_part, part_index)
                )

                for child_node in node.getChildNodes(BODY, part_index, place_in_part=AFTER):
                    envelope.apply_transfer_matrix(
                        self.matrix_factory(child_node, envelope.sync_part)
                    )

            for child_node in node.getChildNodes(EXIT):
                envelope.apply_transfer_matrix(
                    self.matrix_factory(child_node, envelope.sync_part)
                )