import math
import numpy as np

from ..core.bunch import Bunch


class Envelope2D:
    def __init__(self, bunch: Bunch, cov_matrix: np.ndarray = None) -> None:
        self.bunch = bunch

        self.cov_matrix = cov_matrix
        if self.cov_matrix is None:
            self.cov_matrix = np.eye(4)
