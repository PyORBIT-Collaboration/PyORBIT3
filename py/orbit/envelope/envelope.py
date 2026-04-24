import math

import numpy as np


class Envelope:
    def __init__(self, cov_matrix: np.ndarray = None, mean: np.ndarray = None) -> None:
        if mean is None:
            mean = np.zeros(6)

        if cov_matrix is None:
            cov_matrix = np.eye(6)

        self.matrix = np.zeros((7, 7))
        self.matrix[0:6, 0:6] = cov_matrix
        self.matrix[0:6, 6] = mean
        self.matrix[6, 0:6] = mean
        self.matrix[6, 6] = 1.0

    def mean(self) -> np.ndarray:
        return self.matrix[0:6, 6]
    
    def cov(self) -> np.ndarray:
        return self.matrix[0:6, 0:6]
    
    def rms(self) -> np.ndarray:
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
    


    
    