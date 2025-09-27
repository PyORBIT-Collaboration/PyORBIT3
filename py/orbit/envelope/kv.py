"""Envelope model for upright KV distribution."""
import copy
import math
import time
from typing import Callable
from typing import Iterable

# from typing import Self

import numpy as np
import scipy.optimize
from tqdm import tqdm

from orbit.core.bunch import Bunch

from ..lattice import AccActionsContainer
from ..lattice import AccLattice
from ..lattice import AccNode
from ..teapot import TEAPOT_Lattice
from ..teapot import TEAPOT_MATRIX_Lattice
from ..utils import consts

from .nodes import KVEnvelopeTrackerNode
from .lattice_modifications import add_kv_envelope_tracker_nodes
from .utils import calc_cov_twiss
from .utils import bunch_to_numpy
from .utils import get_perveance
from .utils import get_transfer_matrix
from .utils import fit_transfer_matrix
from .utils import build_norm_matrix_from_twiss_2d


class KVEnvelope:
    """Models KV distribution."""
    def __init__(
        self,
        eps_x: float,
        eps_y: float,
        mass: float,
        kin_energy: float,
        line_density: float,
        length: float,
        params: Iterable[float] = None,
    ) -> None:
        """Constructor.

        Args:
            eps_x: RMS emittance in x plane.
            eps_y: RMS emittance in y plane.
            mass: Particle [GeV/c^2].
            kin_energy: Synchronous particle kinetic energy [GeV].
            line_density: Bunch line density [m].
            length: Bunch length [m] (used to convert to Bunch object).
            params: The envelope parameters [cx, cx', cy, cy'].
                The cx and cy parameters represent the envelope extent along the x and y
                axis; cx' and cy' are their derivatives with respect to the distance x.
        """
        self.eps_x = eps_x
        self.eps_y = eps_y
        self.mass = mass
        self.kin_energy = kin_energy

        self.line_density = line_density
        self.length = length
        self.intensity = self.line_density * self.length
        self.perveance = 0.0
        self.set_line_density(line_density)

        self.params = params
        if self.params is None:
            cx = 2.0 * np.sqrt(self.eps_x * 4.0)
            cy = 2.0 * np.sqrt(self.eps_y * 4.0)
            self.params = [cx, 0.0, cy, 0.0]
        self.params = np.array(self.params)

    def set_line_density(self, line_density: float) -> None:
        self.line_density = line_density
        self.intensity = self.line_density * self.length
        self.perveance = get_perveance(self.mass, self.kin_energy, self.line_density)

    def set_length(self, length: float) -> None:
        self.length = length
        self.intensity = self.line_density * self.length

    def set_params(self, params: np.ndarray) -> None:
        self.params = np.copy(params)

    def copy(self):
        return copy.deepcopy(self)

    def cov(self) -> np.ndarray:
        """Return covariance matrix.

        See Table II here: https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.7.024801
        """
        (cx, cxp, cy, cyp) = self.params
        cov_matrix = np.zeros((4, 4))
        cov_matrix[0, 0] = 0.25 * cx**2
        cov_matrix[2, 2] = 0.25 * cy**2
        cov_matrix[1, 1] = 0.25 * cxp**2 + 4.0 * (self.eps_x / cx) ** 2
        cov_matrix[3, 3] = 0.25 * cyp**2 + 4.0 * (self.eps_y / cy) ** 2
        cov_matrix[0, 1] = cov_matrix[1, 0] = 0.25 * cx * cxp
        cov_matrix[2, 3] = cov_matrix[3, 2] = 0.25 * cy * cyp
        return cov_matrix

    def set_cov(self, cov_matrix: np.ndarray) -> None:
        self.eps_x = np.sqrt(np.linalg.det(cov_matrix[0:2, 0:2]))
        self.eps_y = np.sqrt(np.linalg.det(cov_matrix[2:4, 2:4]))
        cx = np.sqrt(4.0 * cov_matrix[0, 0])
        cy = np.sqrt(4.0 * cov_matrix[2, 2])
        cxp = 2.0 * cov_matrix[0, 1] / np.sqrt(cov_matrix[0, 0])
        cyp = 2.0 * cov_matrix[2, 3] / np.sqrt(cov_matrix[2, 2])
        self.set_params([cx, cxp, cy, cyp])

    def twiss(self) -> dict[str, float]:
        cov_matrix = self.cov()
        alpha_x, beta_x, emittance_x = calc_cov_twiss(cov_matrix[0:2, 0:2])
        alpha_y, beta_y, emittance_y = calc_cov_twiss(cov_matrix[2:4, 2:4])

        params = {}
        params["alpha_x"] = alpha_x
        params["alpha_y"] = alpha_y
        params["beta_x"] = beta_x
        params["beta_y"] = beta_y
        params["emittance_x"] = emittance_x
        params["emittance_y"] = emittance_y
        return params

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

    def sample(self, size: int) -> np.ndarray:
        twiss_params = self.twiss()

        alpha_x = twiss_params["alpha_x"]
        alpha_y = twiss_params["alpha_y"]
        beta_x = twiss_params["beta_x"]
        beta_y = twiss_params["beta_y"]
        eps_x = twiss_params["emittance_x"]
        eps_y = twiss_params["emittance_y"]

        V_inv = np.identity(4)
        V_inv[0:2, 0:2] = build_norm_matrix_from_twiss_2d(alpha_x, beta_x)
        V_inv[2:4, 2:4] = build_norm_matrix_from_twiss_2d(alpha_y, beta_y)
        V = np.linalg.inv(V_inv)

        A = np.sqrt(np.diag([eps_x, eps_x, eps_y, eps_y]))

        X = np.random.normal(size=(size, 4))
        X /= np.linalg.norm(X, axis=1)[:, None]
        X /= np.std(X, axis=0)
        X = np.matmul(X, A.T)
        X = np.matmul(X, V.T)
        return X

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

        Args:
            size: Number of macroparticles in the bunch.
                These are the number of  "test"particles not counting the first particle,
                which stores the envelope parameters.
            env: Whether first two particles store envelope parameters.

        Returns:
            Bunch with sampled particles.
        """
        bunch = Bunch()
        bunch.mass(self.mass)
        bunch.getSyncParticle().kinEnergy(self.kin_energy)

        if env:
            (cx, cxp, cy, cyp) = self.params
            bunch.addParticle(cx, cxp, cy, cyp, 0.0, 0.0)

        if size:
            particles = np.zeros((size, 6))
            particles[:, :4] = self.sample(size)
            particles[:, 4] = self.length * np.random.uniform(-0.5, 0.5, size=size)
            for i in range(size):
                bunch.addParticle(*particles[i])

            macrosize = self.intensity / size
            if self.intensity == 0.0:
                macrosize = 1.0
            bunch.macroSize(macrosize)

        return bunch


class KVEnvelopeMonitor:
    def __init__(self, verbose: int = 0) -> None:
        self.verbose = verbose
        self.distance = 0.0
        self._pos_old = 0.0
        self._pos_new = 0.0

        self.history = {}
        for key in [
            "s",
            "xrms",
            "yrms",
        ]:
            self.history[key] = []

    def get_history(self) -> None:
        history = copy.deepcopy(self.history)
        for key in history:
            history[key] = np.array(history[key])
        history["s"] -= history["s"][0]
        return history

    def __call__(self, params_dict: dict) -> None:
        bunch = params_dict["bunch"]
        node = params_dict["node"]

        self._pos_new = params_dict["path_length"]
        if self._pos_old > self._pos_new:
            self._pos_old = 0.0
        self.distance += self._pos_new - self._pos_old
        self._pos_old = self._pos_new

        x_rms = bunch.x(0) * 0.5
        y_rms = bunch.y(0) * 0.5

        self.history["s"].append(self.distance)
        self.history["xrms"].append(x_rms)
        self.history["yrms"].append(y_rms)

        if self.verbose:
            print("s={:0.3f} x_rms={:0.2f}, y_rms={:0.2f}".format(self.distance, x_rms, y_rms))


class KVEnvelopeTracker:
    def __init__(self, lattice: AccLattice, path_length_max: float = None) -> None:
        self.lattice = lattice
        self.nodes = self.add_nodes(
            path_length_min=1.00e-06,
            path_length_max=path_length_max,
        )

        # Lower bounds on envelope parameters
        self.lb = np.zeros(4)
        self.lb[0] = +1.00 - 12
        self.lb[1] = -np.inf
        self.lb[2] = +1.00e-12
        self.lb[3] = -np.inf

        # Upper bounds on envelope parameters
        self.ub = np.zeros(4)
        self.ub[0] = np.inf
        self.ub[1] = np.inf
        self.ub[2] = np.inf
        self.ub[3] = np.inf

    def add_nodes(
        self,
        path_length_min: float,
        path_length_max: float,
    ) -> list[KVEnvelopeTrackerNode]:

        self.nodes = add_kv_envelope_tracker_nodes(
            lattice=self.lattice,
            path_length_min=path_length_min,
            path_length_max=path_length_max,
            perveance=0.0,  # will update based on envelope
            eps_x=1.0,  # will update based on envelope
            eps_y=1.0,  # will update based on envelope
        )
        return self.nodes

    def toggle_nodes(self, setting: bool) -> None:
        for node in self.nodes:
            node.active = setting

    def update_nodes(self, envelope: KVEnvelope) -> None:
        for node in self.nodes:
            node.setPerveance(envelope.perveance)
            node.setEmittances(envelope.eps_x, envelope.eps_y)

    def track(
        self,
        envelope: KVEnvelope,
        periods: int = 1,
        history: bool = False,
    ) -> KVEnvelope:
        self.update_nodes(envelope)

        monitor = KVEnvelopeMonitor()
        action_container = AccActionsContainer()
        if history:
            action_container.addAction(monitor, AccActionsContainer.ENTRANCE)
            action_container.addAction(monitor, AccActionsContainer.EXIT)

        bunch = envelope.to_bunch()        
        for period in range(periods):
            self.lattice.trackBunch(bunch, actionContainer=action_container)

        envelope.from_bunch(bunch)

        if history:
            history = monitor.get_history()
            return (envelope, history)

        return envelope

    def track_particles(self, envelope: KVEnvelope, particles: np.ndarray = None) -> tuple[KVEnvelope, np.ndarray]:
        self.update_nodes(envelope)

        bunch = envelope.to_bunch()
        for i in range(particles.shape[0]):
            bunch.addParticle(*particles[i])
        
        self.lattice.trackBunch(bunch)

        envelope.from_bunch(bunch)

        particles_out = []
        for i in range(1, bunch.getSize()):
            x = bunch.x(i)
            y = bunch.y(i)
            z = bunch.z(i)
            xp = bunch.xp(i)
            yp = bunch.yp(i)
            dE = bunch.dE(i)
            particles_out.append([x, xp, y, yp, z, dE])
        particles_out = np.array(particles_out)

        return (envelope, particles_out)
        
    def transfer_matrix(self, envelope: KVEnvelope) -> np.ndarray:
        bunch = envelope.to_bunch(size=0, env=True)

        if envelope.perveance == 0:
            self.toggle_nodes(False)
            matrix = get_transfer_matrix(self.lattice, bunch)
            self.toggle_nodes(True)
            return matrix

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

        X = bunch_to_numpy(bunch)
        X = X[:, (0, 1, 2, 3)]
        X = X[1:, :]  # ignore first particle, which tracks envelope parameters

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

    def match_zero_sc(self, envelope: KVEnvelope) -> None:
        self.toggle_nodes(False)
        bunch = envelope.to_bunch(size=0, env=False)
        matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice, bunch)
        lattice_params = matrix_lattice.getRingParametersDict()
        self.toggle_nodes(True)

        alpha_x = lattice_params["alpha x"]
        alpha_y = lattice_params["alpha y"]
        beta_x = lattice_params["beta x [m]"]
        beta_y = lattice_params["beta y [m]"]
        envelope.set_twiss(alpha_x, beta_x, alpha_y, beta_y)

    def match(
        self, envelope: KVEnvelope, periods: int = 1, method: str = "least_squares", **kwargs
    ) -> None:
        if envelope.perveance == 0.0:
            return self.match_zero_sc(envelope)

        def loss_function(params: np.ndarray) -> np.ndarray:
            envelope.set_params(params)

            loss = 0.0
            for period in range(periods):
                self.track(envelope)
                residuals = envelope.params - params
                residuals = 1000.0 * residuals
                loss += np.mean(np.abs(residuals))
            return loss / float(periods)

        if method == "least_squares":
            kwargs.setdefault("xtol", 1.00e-12)
            kwargs.setdefault("ftol", 1.00e-12)
            kwargs.setdefault("gtol", 1.00e-12)
            kwargs.setdefault("verbose", 2)

            result = scipy.optimize.least_squares(
                loss_function, envelope.params.copy(), bounds=(self.lb, self.ub), **kwargs
            )
            return result
        elif method == "minimize":
            result = scipy.optimize.minimize(
                loss_function,
                envelope.params.copy(),
                bounds=scipy.optimize.Bounds(self.lb, self.ub),
                **kwargs,
            )
        else:
            raise ValueError
