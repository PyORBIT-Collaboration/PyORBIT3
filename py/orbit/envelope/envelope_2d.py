import math
import numpy as np

from ..core.bunch import SyncParticle


class Envelope2D:
    def __init__(self, sync_part: SyncParticle, cov_matrix: np.ndarray = None) -> None:
        self.sync_part = sync_part

        self.cov_matrix = cov_matrix
        if self.cov_matrix is None:
            self.cov_matrix = np.eye(4)
