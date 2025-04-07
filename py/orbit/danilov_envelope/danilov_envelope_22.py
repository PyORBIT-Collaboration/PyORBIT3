"""Envelope model for {2, 0} Danilov distribution (KV distribution)."""
import copy
import time
from typing import Callable
from typing import Iterable
from typing import Self

import numpy as np
import scipy.optimize
from tqdm import tqdm

from orbit.core.bunch import Bunch

from ..bunch_generators import KVDist2D
from ..bunch_generators import TwissContainer
from ..lattice import AccLattice
from ..lattice import AccNode
from ..teapot import TEAPOT_Lattice
from ..teapot import TEAPOT_MATRIX_Lattice
from ..utils import consts

# from .danilove_envelope_
from .utils import calc_twiss_2d
from .utils import get_bunch_coords
from .utils import get_perveance
from .utils import get_transfer_matrix
from .utils import fit_transfer_matrix


class DanilovEnvelope22:
    """Represents envelope of {2, 2} Danilov distribution (KV distribution with zero emittance in one plane)."""
        
    def set_intensity(self, intensity: int) -> None:
        self.intensity = intensity
        self.line_density = intensity / self.length
        self.perveance = get_perveance(self.mass, self.kin_energy, self.line_density)

    def set_length(self, length: float) -> None:
        self.length = length
        self.set_intensity(self.intensity)

    def set_params(self, params: np.ndarray) -> None:
        self.params = np.copy(params)        
        
    def copy(self) -> Self:
        return copy.deepcopy(self)
        
    def cov(self) -> np.ndarray:
        """Return covariance matrix.
        
        See Table II here: https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.7.024801.
        Note the typo for <xx'> and <yy>, the first time should be rx'^2 not rx^2.
        """
        raise NotImplementedError
    
    def set_cov(self, cov_matrix: np.ndarray) -> None:
        """Set envelope parameters from covariance matrix."""
        raise NotImplementedError
    
    def twiss(self) -> dict[str, float]:
        """Return (alpha_x, beta_x, alpha_y, beta_y)."""
        cov_matrix = self.cov()

        alpha_x, beta_x, emittance_x = calc_twiss_2d(cov_matrix[0:2, 0:2])
        alpha_y, beta_y, emittance_y = calc_twiss_2d(cov_matrix[2:4, 2:4])
        
        results = {}
        results["alpha_x"] = alpha_x
        results["alpha_y"] = alpha_y
        results["beta_x"] = beta_x
        results["beta_y"] = beta_y
        results["emittance_x"] = emittance_x
        results["emittance_y"] = emittance_y
        return results
    
    def set_twiss(self, alpha_x: float = None, beta_x: float = None, alpha_y: float = None, beta_y: float = None) -> None:
        raise NotImplementedError
        
    def from_bunch(self, bunch: Bunch) -> np.ndarray:
        """Set the envelope parameters from a Bunch object."""
        raise NotImplementedError

    def to_bunch(self, size: int = 0, env: bool = True) -> Bunch:
        """Create Bunch object from envelope parameters.
        
        Parameters
        ----------
        size : int
            Number of macroparticles in the bunch. These are the number of "test"
            particles not counting the first particle, which stores the envelope
            parameters.
        env : bool
            If False, do not store the envelope parameters as the first particle.

        Returns
        -------
        Bunch
        """
        raise NotImplementedError


class DanilovEnvelopeTracker22:
    def __init__(self, lattice: AccLattice, path_length_min: float= 1.00e-06) -> None:
        raise NotImplementedError