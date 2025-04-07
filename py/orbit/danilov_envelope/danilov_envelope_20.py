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
from ..bunch_generators import GaussDist2D
from ..bunch_generators import WaterBagDist2D
from ..bunch_generators import TwissContainer
from ..lattice import AccActionsContainer
from ..lattice import AccLattice
from ..lattice import AccNode
from ..teapot import TEAPOT_Lattice
from ..teapot import TEAPOT_MATRIX_Lattice
from ..utils import consts

from .danilov_envelope_solver_nodes import DanilovEnvelopeSolverNode20
from .danilov_envelope_solver_lattice_modifications import add_danilov_envelope_solver_nodes_20
from .utils import calc_twiss_2d
from .utils import get_bunch_coords
from .utils import get_perveance
from .utils import get_transfer_matrix
from .utils import fit_transfer_matrix


class DanilovEnvelope20:
    """Represents envelope of {2, 0} Danilov distribution (KV distribution).

    Attributes
    ----------
    params : ndarray, shape(4,)
        The envelope parameters [cx, cx', cy, cy']. The cx and cy parameters 
        represent the envelope extent along the x and y axis; cx' and cy' are
        their derivatives with respect to the distance x.
    eps_x : float
        The rms emittance of the x-x' distribution: sqrt(<xx><x'x'> - <xx'><xx'>).
    eps_y : float
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
        eps_x: float,
        eps_y: float,
        mass: float,
        kin_energy: float,
        length: float,
        intensity: int,
        params: Iterable[float] = None,
    ) -> None:
        self.eps_x = eps_x
        self.eps_y = eps_y
        self.mass = mass
        self.kin_energy = kin_energy

        self.length = length
        self.line_density = None
        self.perveance = None
        self.set_intensity(intensity)

        self.params = params
        if self.params is None:
            cx = 2.0 * np.sqrt(self.eps_x * 4.0)
            cy = 2.0 * np.sqrt(self.eps_y * 4.0)
            self.params = [cx, 0.0, cy, 0.0]
        self.params = np.array(self.params)
        
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
        
        See Table II here: https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.7.024801
        Note typo for <x'x'>: first term should be rx'^2, not rx^2.
        """
        (cx, cxp, cy, cyp) = self.params
        cov_matrix = np.zeros((4, 4))
        cov_matrix[0, 0] = 0.25 * cx**2
        cov_matrix[2, 2] = 0.25 * cy**2
        cov_matrix[1, 1] = 0.25 * cxp**2 + 4.0 * (self.eps_x / cx)**2
        cov_matrix[3, 3] = 0.25 * cyp**2 + 4.0 * (self.eps_y / cy)**2
        cov_matrix[0, 1] = cov_matrix[1, 0] = 0.25 * cx * cxp
        cov_matrix[2, 3] = cov_matrix[3, 2] = 0.25 * cy * cyp
        return cov_matrix
    
    def set_cov(self, cov_matrix: np.ndarray) -> None:
        """Set envelope parameters from covariance matrix."""
        self.eps_x = np.sqrt(np.linalg.det(cov_matrix[0:2, 0:2]))
        self.eps_y = np.sqrt(np.linalg.det(cov_matrix[2:4, 2:4]))
        cx = np.sqrt(4.0 * cov_matrix[0, 0])
        cy = np.sqrt(4.0 * cov_matrix[2, 2])
        cxp = 2.0 * cov_matrix[0, 1] / np.sqrt(cov_matrix[0, 0])
        cyp = 2.0 * cov_matrix[2, 3] / np.sqrt(cov_matrix[2, 2])
        params = np.array([cx, cxp, cy, cyp])
        self.set_params(params)

    def rms(self) -> np.ndarray:
        return np.sqrt(np.diagonal(self.cov()))
    
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
    
    def set_twiss(
        self, 
        alpha_x: float = None, 
        beta_x: float = None, 
        alpha_y: float = None, 
        beta_y: float = None,
    ) -> None:
        twiss_params = self.twiss()
        if alpha_x is None:
            alpha_x = twiss_params["alpha_x"]
        if alpha_y is None:
            alpha_y = twiss_params["alpha_y"]
        if beta_x is None:
            beta_x = twiss_params["beta_x"]
        if beta_y is None:
            beta_y = twiss_params["beta_y"]
        
        gamma_x = (1.0 + alpha_x**2) / beta_x
        gamma_y = (1.0 + alpha_y**2) / beta_y
        cov_matrix = np.zeros((4, 4))
        cov_matrix[0, 0] = beta_x * self.eps_x
        cov_matrix[2, 2] = beta_y * self.eps_y
        cov_matrix[1, 1] = gamma_x * self.eps_x
        cov_matrix[3, 3] = gamma_y * self.eps_y
        cov_matrix[0, 1] = cov_matrix[1, 0] = -alpha_x * self.eps_x
        cov_matrix[2, 3] = cov_matrix[3, 2] = -alpha_y * self.eps_y
        self.set_cov(cov_matrix)
        
    def sample(self, size: int, dist: str = "kv") -> np.ndarray:
        twiss_params = self.twiss()
        twiss_x = TwissContainer(
            twiss_params["alpha_x"], 
            twiss_params["beta_x"], 
            twiss_params["emittance_x"],
        )
        twiss_y = TwissContainer(
            twiss_params["alpha_y"], 
            twiss_params["beta_y"], 
            twiss_params["emittance_y"],
        )

        if dist == "kv":
            dist = KVDist2D(twiss_x, twiss_y)
        elif dist == "gaussian":
            dist = GaussDist2D(twiss_x, twiss_y)
        elif dist == "waterbag":
            dist = WaterBagDist2D(twiss_x, twiss_y)
        else:
            raise ValueError

        samples = np.zeros((size, 6))
        for i in range(size):
            (x, xp, y, yp) = dist.getCoordinates()
            z = np.random.uniform(-0.5 * self.length, 0.5 * self.length)
            samples[i, :] = [x, xp, y, yp, z, 0.0]
        return samples

    def from_bunch(self, bunch: Bunch) -> np.ndarray:
        """Set envelope parameters from Bunch."""
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
        Bunch
        """
        bunch = Bunch()
        bunch.mass(self.mass)
        bunch.getSyncParticle().kinEnergy(self.kin_energy)

        if env:
            (cx, cxp, cy, cyp) = self.params
            bunch.addParticle(cx, cxp, cy, cyp, 0.0, 0.0)

        if size:
            samples = self.sample(size)
            for i in range(size):
                bunch.addParticle(*samples[i])

            macrosize = self.intensity / size
            if self.intensity == 0.0:
                macrosize = 1.0
            bunch.macroSize(macrosize)

        return bunch
    

class DanilovEnvelopeMonitor20:
    def __init__(self, verbose: int = 0) -> None:
        self.verbose = verbose
        self.position = 0.0

        self.history = {}
        for key in [
            "s",
            "xrms",
            "yrms",
        ]:
            self.history[key] = []
    
    def __call__(self, params_dict: dict) -> None:   
        bunch = params_dict["bunch"]
        node = params_dict["node"]

        if params_dict["path_length"] >= self.position:
            self.position = params_dict["path_length"]
        else:
            self.position = self.position + params_dict["path_length"]

        x_rms = bunch.x(0) * 0.5
        y_rms = bunch.y(0) * 0.5

        self.history["s"].append(self.position)
        self.history["xrms"].append(x_rms)
        self.history["yrms"].append(y_rms)

        if self.verbose:
            print("s={:0.3f} x_rms={:0.2f}, y_rms={:0.2f}".format(position, x_rms, y_rms))
    

class DanilovEnvelopeTracker20:
    def __init__(self, lattice: AccLattice, path_length_max: float = None) -> None:
        self.lattice = lattice
        self.solver_nodes = self.add_solver_nodes(path_length_min=1.00e-06, path_length_max=path_length_max)

        # Lower bounds on envelope parameters
        self.lb = np.zeros(4)
        self.lb[0] = +1.00-12
        self.lb[1] = -np.inf
        self.lb[2] = +1.00e-12
        self.lb[3] = -np.inf

        # Upper bounds on envelope parameters
        self.ub = np.zeros(4)
        self.ub[0] = np.inf
        self.ub[1] = np.inf
        self.ub[2] = np.inf
        self.ub[3] = np.inf

    def add_solver_nodes(
        self, 
        path_length_min: float, 
        path_length_max: float,
    ) -> list[DanilovEnvelopeSolverNode20]:
        
        self.solver_nodes = add_danilov_envelope_solver_nodes_20(
            lattice=self.lattice,
            path_length_min=path_length_min,
            path_length_max=path_length_max,
            perveance=0.0,  # will update based on envelope
            eps_x=1.0,  # will update based on envelope
            eps_y=1.0,  # will update based on envelope
        )
        return self.solver_nodes
    
    def toggle_solver_nodes(self, setting: bool) -> None:
        for node in self.solver_nodes:
            node.active = setting

    def update_solver_node_parameters(self, envelope: DanilovEnvelope20) -> None:
       for solver_node in self.solver_nodes:
            solver_node.set_perveance(envelope.perveance)
            solver_node.set_emittances(envelope.eps_x * 4.0, envelope.eps_y * 4.0)

    def track(self, envelope: DanilovEnvelope20, periods: int = 1, history: bool = False) -> None | dict[str, np.ndarray]:
        self.update_solver_node_parameters(envelope)

        monitor = DanilovEnvelopeMonitor20()
        action_container = AccActionsContainer()
        if history:
            action_container.addAction(monitor, AccActionsContainer.ENTRANCE)
            action_container.addAction(monitor, AccActionsContainer.EXIT)

        bunch = envelope.to_bunch()
        for period in range(periods):
            self.lattice.trackBunch(bunch, actionContainer=action_container)

        envelope.from_bunch(bunch)

        if history:
            for key in  monitor.history:
                monitor.history[key] = np.array(monitor.history[key])
            return monitor.history
    
    def get_transfer_matrix(self, envelope: DanilovEnvelope20) -> np.ndarray:
        bunch = envelope.to_bunch(size=0, env=False)

        if self.perveance == 0:
            return get_transfer_matrix(self.lattice, bunch)

        step_arr = np.ones(6) * 1.00e-06
        step_reduce = 20.0

        bunch.addParticle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(step_arr[0] / step_reduce, 0.0, 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, step_arr[1] / step_reduce, 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, 0.0, step_arr[2] / step_reduce, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, 0.0, 0.0, step_arr[3] / step_reduce, 0.0, 0.0)
        bunch.addParticle(step_arr[0], 0.0, 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, step_arr[1], 0.0, 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, 0.0, step_arr[2], 0.0, 0.0, 0.0)
        bunch.addParticle(0.0, 0.0, 0.0, step_arr[3], 0.0, 0.0)

        self.lattice.trackBunch(bunch)
        
        X = get_bunch_coords(bunch)
        X = X[:, (0, 1, 2, 3)]
        X = X[1:, :]

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
    
    def match_zero_sc(self, envelope: DanilovEnvelope20) -> None:
        self.toggle_solver_nodes(False)
        bunch = envelope.to_bunch(size=0, env=False)
        matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice, bunch)
        lattice_params = matrix_lattice.getRingParametersDict()
        self.toggle_solver_nodes(True)
        
        alpha_x = lattice_params["alpha x"]
        alpha_y = lattice_params["alpha y"]
        beta_x = lattice_params["beta x [m]"]
        beta_y = lattice_params["beta y [m]"]
        envelope.set_twiss(alpha_x, beta_x, alpha_y, beta_y)
                        
    def match(self, envelope: DanilovEnvelope20, method: str = "least_squares", **kwargs) -> None:
        if envelope.perveance == 0.0:
            return self.match_zero_sc(envelope)
                
        def loss_function(params: np.ndarray) -> np.ndarray:   
            envelope.set_params(params)
            self.track(envelope)
            residuals = envelope.params - params
            residuals = 1000.0 * residuals        
            loss = np.mean(np.abs(residuals))
            return loss

        if method == "least_squares":
            kwargs.setdefault("xtol", 1.00e-12)
            kwargs.setdefault("ftol", 1.00e-12)
            kwargs.setdefault("gtol", 1.00e-12)
            kwargs.setdefault("verbose", 2)

            result = scipy.optimize.least_squares(
                loss_function,
                envelope.params.copy(),
                bounds=(self.lb, self.ub),
                **kwargs
            )
            return result
        else:
            raise ValueError