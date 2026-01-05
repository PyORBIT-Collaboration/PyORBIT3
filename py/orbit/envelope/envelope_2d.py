import math
import numpy as np

from ..core.bunch import Bunch


class Envelope2D:
    def __init__(self, bunch: Bunch, cov_matrix: np.ndarray = None) -> None:
        self.bunch = bunch

        self.cov_matrix = cov_matrix
        if self.cov_matrix is None:
            self.cov_matrix = np.eye(4)

    def track(self, matrix: np.ndarray) -> None:
        """Evolve covariance matrix: S -> M S M^T."""
        S = self.cov_matrix
        M = matrix[:4, :4]
        self.cov_matrix = np.linalg.multi_dot([M, S, M.T])

    def spaceChargeMatrix(self) -> np.ndarray:
        """Return transfer matrix for linear space charge kick."""
        raise NotImplementedError()