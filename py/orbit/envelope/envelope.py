import math

import numpy as np

from ..lattice import AccNode
from ..lattice import AccLattice

from .matrix import MatrixFactory


ENTRANCE = AccNode.ENTRANCE
BODY = AccNode.BODY
EXIT = AccNode.EXIT

BEFORE = AccNode.BEFORE
AFTER = AccNode.AFTER


class Envelope:
    """Represents beam envelope and centroid in 6D phase space.
    
    Attributes:
        matrix: 7x7 matrix containing 6x6 covariance matrix and 6x1 centroid vector.
    """
    def __init__(self, cov_matrix: np.ndarray = None, centroid: np.ndarray = None) -> None:
        """Constructor.
        
        Args:
            cov_matrix: 6x6 covariance matrix.
            centroid: 6x1 centroid vector.
        """
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
        """Return centroid vector."""
        return self.matrix[0:6, 6]
    
    def cov(self) -> np.ndarray:
        """Return covariance matrix."""
        return self.matrix[0:6, 0:6]
    
    def rms(self) -> np.ndarray:
        """Return rms beam sizes."""
        return np.sqrt(np.diag(self.cov()))

    def apply_transfer_matrix(self, transfer_matrix: np.ndarray) -> None:
        """Linear propagation of beam matrix.
        
        Args:
            transfer_matrix: 7x7 transfer matrix.
        """
        self.matrix = np.linalg.multi_dot([transfer_matrix, self.matrix, transfer_matrix.T])

    def space_charge_matrix(self, length: float) -> np.ndarray:
        """Return transfer matrix from linear space charge kick."""
        raise NotImplementedError()
    

def get_matrix(node: AccNode) -> np.ndarray:
    """Return transfer matrix for given node."""

