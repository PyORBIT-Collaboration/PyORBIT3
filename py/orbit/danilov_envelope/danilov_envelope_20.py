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

from .utils import calc_twisss_2d
from .utils import get_transfer_matrix


class DanilovEnvelope20:
    """Represents envelope of {2, 0} Danilov distribution (KV distribution).

    Attributes
    ----------
    params : ndarray, shape(4,)
        The envelope parameters [cx, cx', cy, cy']. The cx and cy parameters 
        represent the envelope extent along the x and y axis; cx' and cy' are
        their derivatives with respect to the distance x.
    eps_x_rms : float
        The rms emittance of the x-x' distribution: sqrt(<xx><x'x'> - <xx'><xx'>).
    eps_y_rms : float
        The rms emittance of the y-y' distribution: sqrt(<yy><y'y'> - <yy'><yy'>).
    mass : float
        Particle [GeV/c^2].
    kin_energy : float
        Particle kinetic energy [GeV].
    intensity : float
        Bunch intensity (number of particles).
    length : float
        Bunch length [m].
    perveance : float
        Dimensionless beam perveance.
    """
    def __init__(
        self,
        eps_x_rms: float,
        eps_y_rms: float,
        mass: float,
        kin_energy: float,
        length: float, 1.0,
        intensity: int,
        params: Iterable[float] = None,
    ) -> None:
        self.eps_x_rms = eps_x_rms
        self.eps_y_rms = eps_y_rms
        self.eps_x = None
        self.eps_y = None
        self.set_rms_emittances(eps_x_rms, eps_y_rms)

        self.mass = mass
        self.kin_energy = kin_energy
        self.length = length

        self.line_density = None
        self.perveance = None
        self.set_intensity(intensity)

        self.params = params
        if self.params is None:
            cx = 2.0 * np.sqrt(self.eps_x)
            cy = 2.0 * np.sqrt(self.eps_y)
            self.params = [cx, 0.0, cy, 0.0]

        self.params = np.array(self.params)
        
    def set_rms_emittances(self, eps_x_rms: float, eps_y_rms: float) -> None:
        self.eps_x_rms = eps_x_rms
        self.eps_y_rms = eps_y_rms
        self.eps_x = 4.0 * eps_x_rms
        self.eps_y = 4.0 * eps_y_rms
        
    def set_intensity(self, intensity: int) -> None:
        self.intensity = intensity
        self.line_density = intensity / self.length
        self.perveance = utils.get_perveance(self.mass, self.kin_energy, self.line_density)

    def set_length(self, length: float) -> None:
        self.length = length
        self.set_intensity(self.intensity)

    def set_params(self, params):
        self.params = params        
        
    def copy(self) -> Self:
        return copy.deepcopy(self)
        
    def cov(self) -> np.ndarray:
        """Return covariance matrix.
        
        See Table II here: https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.7.024801.
        Note the typo for <xx'> and <yy>, the first time should be rx'^2 not rx^2.
        """
        (cx, cxp, cy, cyp) = self.params
        Sigma = np.zeros((4, 4))
        Sigma[0, 0] = cx ** 2
        Sigma[2, 2] = cy ** 2
        Sigma[1, 1] = cxp**2 + (self.eps_x / cx)**2
        Sigma[3, 3] = cyp**2 + (self.eps_y / cy)**2
        Sigma[0, 1] = Sigma[1, 0] = cx * cxp
        Sigma[2, 3] = Sigma[3, 2] = cy * cyp
        Sigma = Sigma * 0.25
        return Sigma
    
    def set_cov(self, cov_matrix: np.ndarray) -> None:
        """Set envelope parameters from covariance matrix."""
        eps_x = np.sqrt(np.linalg.det(cov_matrix[0:2, 0:2]))
        eps_y = np.sqrt(np.linalg.det(cov_matrix[2:4, 2:4]))
        cx = 2.0 * np.sqrt(cov_matrix[0, 0])
        cy = 2.0 * np.sqrt(cov_matrix[2, 2])
        cxp = 4.0 * cov_matrix[0, 1] / cx
        cyp = 4.0 * cov_matrix[2, 3] / cy
        self.set_params([cx, cxp, cy, cyp])
        self.set_rms_emittances(eps_x, eps_y)
    
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
        
    def from_bunch(self, bunch: Bunch) -> np.ndarray:
        """Set the envelope parameters from a Bunch object."""
        self.params = np.zeros(4)
        self.params[0] = bunch.x(0)
        self.params[1] = bunch.xp(0)
        self.params[2] = bunch.y(0)
        self.params[3] = bunch.yp(0)
        return self.params

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
        bunch: Bunch object
            The bunch representing the distribution of size 2 + n_parts
            (unless `no_env` is True).
        params_dict : dict
            The dictionary of parameters for the bunch.
        """
        bunch = Bunch()
        bunch.mass(self.mass)
        bunch.getSyncParticle().kinEnergy(self.kin_energy)

        params_dict = {"bunch": bunch}
        if env:
            cx, cxp, cy, cyp = self.params
            bunch.addParticle(cx, cxp, cy, cyp, 0.0, 0.0)

        twiss_params = self.twiss()
        dist = KVDist2D(
            TwissContainer(twiss_params["alpha_x"], twiss_params["beta_x"], twiss_params["emittance_x"]),
            TwissContainer(twiss_params["alpha_y"], twiss_params["beta_y"], twiss_params["emittance_y"]),
        )
        if size:
            for i in range(size):
                (x, xp, y, yp) = dist.getCoordinates()
                z = np.random.uniform(0.0, self.length)
                bunch.addParticle(x, xp, y, yp, z, 0.0)

            macrosize = self.intensity / size
            if self.intensity == 0.0:
                macrosize = 1.0
            bunch.macroSize(macrosize)

        return bunch

    def track(self, lattice: AccLattice, periods: int = 1) -> None:
        """Track the envelope through the lattice.
        
        If `size > 0` test particles will be tracked which receive linear space charge
        kicks based on the envelope parameters. These particles are sampled from a KV
        distribution.
        """
        bunch = self.to_bunch(size)
        for _ in range(periods):
            lattice.trackBunch(bunch)
        self.from_bunch(bunch)

    def track_store_params(self, lattice: AccLattice, periods: int = 1) -> None:
        """Track and return the turn-by-turn envelope parameters."""
        params_tbt = [self.params]
        for _ in range(periods):
            self.track(lattice)
            params_tbt.append(self.params)
        return params_tbt

    def get_transfer_matrix(self, lattice: AccLattice) -> np.ndarray:
        """Compute effective transfer matrix with linear space charge included.

        The method is taken from /src/teapot/MatrixGenerator.cc.
        
        Parameters
        ----------
        lattice : TEAPOT_Lattice object
            The lattice may have envelope solver nodes. These nodes should
            track the beam envelope using the first two particles in the bunch,
            then use these to apply the appropriate linear space charge kicks
            to the rest of the particles.

        Returns
        -------
        M : ndarray, shape (4, 4)
            Effective transfer matrix with linear space charge included.
        """
        # Use TEAPOT_MATRIX_Lattice if perveance is zero.
        if self.perveance == 0:
            M = get_transfer_matrix(lattice, mass=self.mass, kin_energy=self.kin_energy)
            return M[:4, :4]

        # The envelope parameters will change if the beam is not matched to
        # the lattice, so make a copy of the intial state.
        env = self.copy()

        step_arr_init = np.full(6, 1.00e-06)
        step_arr = np.copy(step_arr_init)
        step_reduce = 20.0
 
        bunch = env.to_bunch()
        bunch.addParticle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(step_arr[0] / step_reduce, 0.0, 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, step_arr[1] / step_reduce, 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, 0.0, step_arr[2] / step_reduce, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, 0.0, 0.0, step_arr[3] / step_reduce, 0.0, 0.0)
        bunch.addParticle(step_arr[0], 0.0, 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, step_arr[1], 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, 0.0, step_arr[2], 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, 0.0, 0.0, step_arr[3], 0.0, 0.0)

        lattice.trackBunch(bunch, params_dict)
        X = [
            [bunch.x(i), bunch.xp(i), bunch.y(i), bunch.yp(i)]
            for i in range(1, bunch.getSize())
        ]
        X = np.array(X)            
        M = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                x1 = step_arr[i] / step_reduce
                x2 = step_arr[i]
                y0 = X[0, j]
                y1 = X[i + 1, j]
                y2 = X[i + 1 + 4, j]
                M[j, i] = ((y1 - y0) * x2 * x2 - (y2 - y0) * x1 * x1) / (x1 * x2 * (x2 - x1))
        return M

    def match_bare(self, lattice, solver_nodes=None):
        """Match to the lattice without space charge.

        Parameters
        ----------
        lattice : TEAPOT_Lattice object
            The lattice in which to match. If envelope solver nodes nodes are
            in the lattice, a list of these nodes needs to be passed as the
            `solver_nodes` parameter so they can be turned off/on.
        solver_nodes : list, optional
            List of nodes which are sublasses of SC_Base_AccNode. If provided,
            all space charge nodes are turned off, then the envelope is matched
            to the bare lattice, then all space charge nodes are turned on.

        Returns
        -------
        ndarray, shape (4,)
            The matched envelope parameters.
        """
        if solver_nodes is not None:
            for node in solver_nodes:
                node.active = False

        M = self.get_transfer_matrix(lattice)
        tmat = TransferMatrix(M)
        if not tmat.stable:
            print("WARNING: transfer matrix is not stable.")

        tmat = TransferMatrixCourantSnyder(M)
        tmat.analyze()
        alpha_x = tmat.params["alpha_x"]
        alpha_y = tmat.params["alpha_y"]
        beta_x = tmat.params["beta_x"]
        beta_y = tmat.params["beta_y"]
        sig_xx = self.eps_x_rms * beta_x 
        sig_yy = self.eps_y_rms * beta_y
        sig_xxp = -self.eps_x_rms * alpha_x
        sig_yyp = -self.eps_y_rms * alpha_y
        cx = 2.0 * np.sqrt(sig_xx)
        cy = 2.0 * np.sqrt(sig_yy)
        cxp = 4.0 * sig_xxp / cx
        cyp = 4.0 * sig_yyp / cy
        self.set_params([cx, cxp, cy, cyp])

        if solver_nodes is not None:
            for node in solver_nodes:
                node.active = True
                
        return self.params
        
    def match_lsq(self, lattice, **kws):
        """Compute matched envelope using scipy.least_squares optimizer.

        Parameters
        ----------
        lattice : TEAPOT_Lattice object
            The lattice to match into. The solver nodes should already be in place.
        **kws
            Key word arguments passed to `scipy.optimize.least_squares`.
            
        Returns
        -------
        result : scipy.optimize.OptimizeResult
            See scipy documentation.
        """        
        def cost_function(x):   
            x0 = x.copy()
            self.params = x
            self.track(lattice)
            residuals = self.params - x0
            residuals = 1000.0 * residuals
            return residuals
        
        kws.setdefault("xtol", 1.00e-12)
        kws.setdefault("ftol", 1.00e-12)
        kws.setdefault("gtol", 1.00e-12)

        lb = [1.00e-12, -np.inf, 1.00e-12, -np.inf]
        result = scipy.optimize.least_squares(
            cost_function,
            self.params.copy(),
            bounds=(lb, np.inf),
            **kws
        )
        self.params = result.x
        return result
    
    def match_lsq_ramp_intensity(self, lattice, solver_nodes=None, n_steps=10, **kws):
        if self.perveance == 0.0:
            return
        verbose = kws.get("verbose", 0)
        max_intensity = self.intensity
        for intensity in np.linspace(0.0, max_intensity, n_steps):
            self.set_intensity(intensity)
            for solver_node in solver_nodes:
                solver_node.set_perveance(self.perveance)
            result = self.match_lsq(lattice, **kws)
            self.set_params(result.x)
            if verbose > 0:
                print("intensity = {:.2e}".format(intensity))
        return result