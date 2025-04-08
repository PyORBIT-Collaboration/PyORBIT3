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

from .utils import get_bunch_coords
from .utils import get_perveance
from .utils import get_transfer_matrix
from .utils import fit_transfer_matrix
from .utils import rotation_matrix_xy

from .transfer_matrix import norm_matrix_from_twiss_cs
from .transfer_matrix import norm_matrix_from_twiss_lb_one_mode


def cov_matrix_to_vector(cov_matrix: np.ndarray) -> np.ndarray:
    """Return array of 10 unique elements of covariance matrix."""
    return cov_matrix[np.triu_indices(4)]


def env_matrix_to_vector(env_matrix: np.ndarray) -> np.ndarray:
    """Return list of envelope parameters from envelope matrix."""
    return env_matrix.ravel()


def env_vector_to_matrix(env_vector: np.ndarray) -> np.ndarray:
    return env_vector.reshape(4, 2)


def env_params_from_bunch(bunch: Bunch) -> np.ndarray:
    a = bunch.x(0)
    b = bunch.x(1)
    e = bunch.y(0)
    f = bunch.y(1)
    ap = bunch.xp(0)
    bp = bunch.xp(1)
    ep = bunch.yp(0)
    fp = bunch.yp(1)
    params = [a, b, ap, bp, e, f, ep, fp]
    params = np.array(params)
    return params


class DanilovEnvelope22:
    """Represents envelope of {2, 2} Danilov distribution (KV distribution with zero emittance in one plane).

    Attributes
    ----------
    params: np.ndarray, shape (8,)
        The envelope parameters [a, b, a', b', e, f, e', f']. The coordinates
        of a particle on the beam envelope are parameterized as
            x = a*cos(psi) + b*sin(psi), x' = a'*cos(psi) + b'*sin(psi),
            y = e*cos(psi) + f*sin(psi), y' = e'*cos(psi) + f'*sin(psi),
        where 0 <= psi <= 2pi.
    intrinsic_emittance: float
        Nonzero rms emittance [m * rad].
    mode : int 
        Mode number corresponding to nonzero emittance. In certain coordinates this determines
        the sign of the angular momentum.
    eps_x_frac : float
        Determines horizontal emittance relative to intrinsic emittance:
        eps_x = eps_x_frac * intrinsic_emittance. This is determined by the phase
        of the mode eigenvector.
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
        intrinsic_emittance: float,
        eps_x_frac: float,
        mass: float,
        kin_energy: float,
        length: float,
        intensity: float = 0.0,
        mode: int = 1,
        params: np.ndarray = None,
    ) -> None:
        self.intrinsic_emittance = intrinsic_emittance
        self.mode = mode
        self.eps_x_frac = eps_x_frac
        self.mass = mass
        self.kin_energy = kin_energy

        self.length = length
        self.line_density = None
        self.perveance = None
        self.set_intensity(intensity)

        if params is None:
            eps_x = self.intrinsic_emittance * self.eps_x_frac
            eps_y = self.intrinsic_emittance * (1.0 - self.eps_x_frac)
            rx = 2.0 * np.sqrt(eps_x)
            ry = 2.0 * np.sqrt(eps_y)
            if mode == 1:
                params = np.array([rx, 0, 0, rx, 0, -ry, +ry, 0])
            else:
                params = np.array([rx, 0, 0, rx, 0, +ry, -ry, 0])

        self.params = params
        self.set_params(self.params)

    def copy(self) -> Self:
        return copy.deepcopy(self)

    def set_intensity(self, intensity: int) -> None:
        self.intensity = intensity
        self.line_density = intensity / self.length
        self.perveance = get_perveance(self.mass, self.kin_energy, self.line_density)

    def set_length(self, length: float) -> None:
        self.length = length
        self.set_intensity(self.intensity)

    def set_params(self, params: np.ndarray) -> None:
        self.params = np.copy(params)    

        eps_x, eps_y = self.projected_emittances()
        self.intrinsic_emittance = eps_x + eps_y
        self.eps_x_frac = eps_x / self.intrinsic_emittance    

    def param_matrix(self) -> np.ndarray:
        """Build envelope matrix.

        The envelope matrix P defined by x = Pc, where x = [x, x', y, y']^T and 
        c = [cos(psi), sin(psi)], with 0 <= psi <= 2pi. This is useful because 
        any transformation to the phase space coordinate vector x is also done to
        P. For example, if x -> Mx, then P -> MP.
        """
        return env_vector_to_matrix(self.params)

    def param_vector(self, axis: int) -> np.ndarray:
        """Return vector of envelope parameters [a, b, a', b', e, f, e', f']."""
        if axis is None:
            return self.params
        param_matrix = self.param_matrix()
        return param_matrix[axis, :]
    
    def transform(self, matrix: np.ndarray) -> None:
        """Linearly transform phase space coordinates."""
        param_matrix = self.param_matrix()
        param_matrix = np.matmul(matrix, param_matrix)
        self.params = np.ravel(param_matrix)

    def rotate(self, angle: float) -> None:
        """Apply clockwise rotation by `angle`` radians in x-y plane."""
        self.transform(rotation_matrix_xy(angle))

    def normalization_matrix(self, kind: str = "2d") -> np.ndarray:
        if kind == "2d":
            norm_matrix = CS.normalization_matrix_from_twiss(*self.twiss_2d())
        elif kind == "4D":
            alpha_lx, beta_lx, alpha_ly, beta_ly, u, nu = self.twiss_4d()
            norm_matrix = LB.normalization_matrix_from_twiss_one_mode(
                alpha_lx=alpha_lx,
                beta_lx=beta_lx,
                alpha_ly=alpha_ly,
                beta_ly=beta_ly,
                u=u,
                nu=nu,
                mode=self.mode,
            )
        else:
            raise ValueError(f"Invalid normalization {kind}")
        return norm_matrix
    
    def normalize(self, kind: str = "2d", scale: bool = False) -> None:
        """Normalize the phase space coordinates.

        Parameters
        ----------
        kind : {"2d", "4d"}
            - "2d": The x-x' and y-y' ellipses will be circles of radius sqrt(eps_x)
                    and sqrt(eps_y), where eps_x and eps_y are the rms projected
                    emittances.
            - "4d": The 4 x 4 covariance matrix becomes diagonal. The x-x' and y-y'
                    ellipses wil be circles of radius radius sqrt(eps_1) and
                    sqrt(eps_2), where eps_1, and eps_2 are the rms intrinsic
                    emittances.
        scale : bool
            Whether to divide by the emittances to scale all coordinates to unit
            variance. This converts all circles to unit circles.
        """
        if kind == "2d":
            self.transform(self.normalization_matrix(kind=kind))
            if scale:
                eps_x, eps_y = self.projected_emittances()
                if eps_x > 0.0:
                    self.params[0:4] /= np.sqrt(4.0 * eps_x)
                if eps_y > 0.0:
                    self.params[4:0] /= np.sqrt(4.0 * eps_y)
            return self.params
        
        elif kind == "4d":
            r_n = np.sqrt(4.0 * self.intrinsic_emittance)
            if self.mode == 1:
                self.params = np.array([r_n, 0, 0, r_n, 0, 0, 0, 0])
            else:
                self.params = np.array([0, 0, 0, 0, 0, r_n, r_n, 0])
            if scale:
                self.params = self.params / r_n
            return self.params
        
        else:
            raise ValueError(f"Invalid normalization {kind}")
        
    def projected_tilt_angle(self, axis: tuple[int, int] = (0, 2)) -> float:
        """Return angle of projected ellipse."""
        a, b = self.param_vector(axis[0])
        e, f = self.param_vector(axis[1])
        return 0.5 * np.arctan2(2 * (a * e + b * f), a**2 + b**2 - e**2 - f**2)
    
    def projected_radii(self, axis: tuple[int, int] = (0, 2)) -> tuple[float, float]:
        """Return semi-major and semi-minor axes of projected ellipse."""
        a, b = self.param_vector(axis[0])
        e, f = self.param_vector(axis[1])
        phi = self.projected_tilt_angle(axis)
        _sin = np.sin(phi)
        _cos = np.cos(phi)
        _sin2 = _sin ** 2
        _cos2 = _cos ** 2
        xx = a**2 + b**2
        yy = e**2 + f**2
        xy = a * e + b * f
        cx = np.sqrt(abs(xx * _cos2 + yy * _sin2 - 2 * xy * _sin * _cos))
        cy = np.sqrt(abs(xx * _sin2 + yy * _cos2 + 2 * xy * _sin * _cos))
        return (cx, cy)
    
    def projected_area(self, axis: tuple[int, int] = (0, 2)) -> float:
        """Return area of projected ellipse."""
        a, b = self.param_vector(axis[0])
        e, f = self.param_vector(axis[1])
        return np.pi * abs(a * f - b * e)

    def cov(self) -> np.ndarray:
        """Return covariance matrix."""
        param_matrix = self.param_matrix()
        return 0.25 * np.matmul(param_matrix, param_matrix.T)
    
    def corr(self) -> np.ndarray:
        """Return correlation matrix."""
        S = self.cov()
        D = np.sqrt(np.diag(S.diagonal()))
        D_inv = np.linalg.inv(D)
        return np.linalg.multi_dot([D_inv, S, D_inv])
    
    def projected_emittances(self) -> tuple[float, float]:
        """Return rms apparent emittances eps_x, eps_y [m * rad]."""
        eps_x = eps_y = 0.0
        cov_matrix = self.cov()
        determinant = np.linalg.det(cov_matrix[0:2, 0:2])
        if determinant > 0.0:
            eps_x = np.sqrt(determinant)
        determinant = np.linalg.det(cov_matrix[2:4, 2:4])
        if determinant > 0.0:
            eps_y = np.sqrt(determinant)
        return (eps_x, eps_y)

    def intrinsic_emittances(self) -> tuple[float, float]:
        """Return rms intrinsic emittances eps1, eps2 [m * rad]."""
        eps_1 = eps_2 = 0.0
        if self.mode == 1:
            eps_1 = self.intrinsic_emittance
        else:
            eps_2 = self.intrinsic_emittance
        return (eps_1, eps_2)
    
    def phases_xy(self) -> tuple[float, float]:
        envelope = self.copy()
        envelope.normalize("2d")
        a, b, ap, bp, e, f, ep, fp = envelope.params
        phase_x = -np.arctan2(ap, a)
        phase_y = -np.arctan2(ep, e)
        if phase_x < 0.0:
            phase_x += 2.0 * np.pi
        if phase_y < 0.0:
            phase_y += 2.0 * np.pi
        return (phase_x, phase_y)
        
    def nu(self) -> float:
        """Return the x-y phase difference (nu) of all particles in the beam.

        The value returned is in the range [0, pi]. This can also be found from
        the equation cos(nu) = r, where r is the x-y correlation coefficient.
        """
        phase_x, phase_y = self.phases_xy()
        nu = abs(phase_y - phase_x)
        if nu < np.pi:
            return nu
        else:
            return 2.0 * np.pi - nu
        
    def u(self) -> float:
        cov_matrix = self.cov
        if self.mode == 1:
            determinant = np.linalg.det(cov_matrix[2:4, 2:4])
            eps_y = 0.0
            if determinant > 0.0:
                eps_y = np.sqrt(determinant)
            u = eps_y / self.intrinsic_emittance
        else:
            determinant = np.linalg.det(cov_matrix[0:2, 0:2])
            eps_x = 0.0
            if determinant > 0.0:
                eps_x = np.sqrt(determinant)
            u = eps_x / self.intrinsic_emittance
        return u
    
    def twiss_2d(self) -> dict[str, float]:
        """Return Twiss parameters in x-x', y-y' planes."""
        eps_x = 0.0
        eps_y = 0.0
        beta_x = np.inf
        beta_y = np.inf
        alpha_x = np.inf
        alpha_y = np.inf

        cov_matrix = self.cov()
        determinant = np.linalg.det(cov_matrix[:2, :2])
        if determinant > 0.0:
            eps_x = np.sqrt(determinant)
            beta_x = cov_matrix[0, 0] / eps_x
            alpha_x = -cov_matrix[0, 1] / eps_x

        determinant = np.linalg.det(cov_matrix[2:, 2:])
        if determinant > 0.0:
            eps_y = np.sqrt(determinant)
            beta_y = cov_matrix[2, 2] / eps_y
            alpha_y = -cov_matrix[2, 3] / eps_y

        return {
            "alpha_x": alpha_x,
            "alpha_y": alpha_y,
            "beta_x": beta_x,
            "beta_y": beta_y,
        }
    
    def twiss_4d(self) -> dict[str, float]:
        """Return Bogacz/Lebedev 4D Twiss parameters for occupied mode."""
        cov_matrix = self.cov()
        beta_lx = cov_matrix[0, 0] / self.intrinsic_emittance
        beta_ly = cov_matrix[2, 2] / self.intrinsic_emittance
        alpha_lx = -cov_matrix[0, 1] / self.intrinsic_emittance
        alpha_ly = -cov_matrix[2, 3] / self.intrinsic_emittance
        u = self.u()
        nu = self.nu()
        return {
            "alpha_lx": alpha_lx,
            "alpha_ly": alpha_ly,
            "beta_lx": beta_lx,
            "beta_ly": beta_ly,
            "u": u,
            "nu": nu,
        }
    
    def set_twiss_2d(
        self, 
        alpha_x: float = None, 
        beta_x: float = None,
        alpha_y: float = None,
        beta_y: float = None,
        eps_x_frac: float = None,
    ) -> None:
        
        twiss_params = self.twiss_2d()
        if alpha_x is None:
            alpha_x = twiss_params["alpha_x"]
        if alpha_y is None:
            alpha_y = twiss_params["alpha_y"]
        if beta_x is None:
            beta_x = twiss_params["beta_x"]  
        if beta_y is None:      
            beta_y = twiss_params["beta_y"]

        V_inv = norm_matrix_from_twiss_cs(alpha_x=alpha_x, beta_x=beta_x, alpha_y=alpha_y, beta_y=beta_y)
        V = np.linalg.inv(V_inv)

        if eps_x_frac is None:
            eps_x_frac = self.eps_x_frac

        eps_x = self.intrinsic_emittance * eps_x_frac
        eps_y = self.intrinsic_emittance * (1.0 - eps_x_frac)
        A = np.sqrt(4.0 * np.diag([eps_x, eps_x, eps_y, eps_y]))
        V = np.matmul(V, A)
        
        self.normalize(method="2d", scale=True)
        self.transform(V)

    def set_twiss_4d(
        self, 
        alpha_lx: float = None,
        beta_lx: float = None, 
        alpha_ly: float = None, 
        beta_ly: float = None, 
        u: float = None, 
        nu: float = None
    ) -> None:
        """Set 4D Twiss parameters.

        Parameters
        ----------
        - alpha_lx : Horizontal alpha function -<xx'> / eps_l.
        - beta_lx : Horizontal beta function <xx> / eps_l.
        - alpha_ly : Vertical alpha_function -<yy'> / eps_l.
        - beta_ly : Vertical beta function -<yy> / eps_l.
        - u : Coupling parameter in range [0, 1] (eps_y / epsl if mode=1, eps_x / epsl if mode=2).
        - nu : The x-y phase difference in range [0, pi].
        """
        twiss_params = self.twiss_4d()
        for i, value in enumerate([alpha_lx, beta_lx, alpha_ly, beta_ly, u, nu]):
            if value is not None:
                twiss_params[i] = value

        V_inv = norm_matrix_from_twiss_lb_one_mode(
            alpha_lx=alpha_lx,
            beta_lx=beta_lx,
            alpha_ly=alpha_ly,
            beta_ly=beta_ly,
            u=u,
            nu=nu,
            mode=1,
        )
        V = np.linalg.inv(V_inv)

        self.normalize(method="4d")
        self.transform(V)

    def set_cov(self, cov_matrix: np.ndarray, verbose: int = 0) -> scipy.optimize.OptimizeResult:
        """Set envelope parameters from covariance matrix."""

        def loss_function(params: np.ndarray, cov_matrix: np.ndarray) -> np.ndarray:
            self.set_params(params)
            loss = cov_matrix_to_vector(cov_matrix - self.cov())
            loss = loss * 1.00e+06
            return loss

        result = scipy.optimize.least_squares(
            loss_function, self.params, args=(cov_matrix,), xtol=1.00e-12, verbose=verbose
        )
        return result
    
    def get_coordinates(self, psi: float = 0.0) -> np.ndarray:
        """Return phase space coordinates of particle on envelope.

        psi is in the range [0, 2 * pi].
        """
        param_matrix = self.param_matrix()
        return np.matmul(param_matrix, [np.cos(psi), np.sin(psi)])
        
    def from_bunch(self, bunch: Bunch) -> None:
        """Set the envelope parameters from a Bunch object."""
        self.set_params(env_params_from_bunch(bunch))

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
            a, b, ap, bp, e, f, ep, fp = self.params
            bunch.addParticle(a, ap, e, ep, 0.0, 0.0)
            bunch.addParticle(b, bp, f, fp, 0.0, 0.0)

        if size:
            samples = self.sample(size)
            for i in range(size):
                bunch.addParticle(*samples[i])

            macrosize = self.intensity / size
            if self.intensity == 0.0:
                macrosize = 1.0
            bunch.macroSize(macrosize)

        return bunch
        
    def sample(self, size: int, dist: str = "kv") -> np.ndarray:
        if dist == "kv":
            psis = np.linspace(0.0, 2.0 * np.pi, size)
            samples = [self.get_coordinates(psi) for psi in psis]
            samples = np.array(samples)
            radii = np.random.uniform(0.0, 1.0, size=size) ** -0.5
            samples = samples * radii[:, None]
            return samples
        elif dist == "gaussian":
            return np.random.multivariate_normal(np.zeros(4), self.cov(), size=size)
        else:
            raise ValueError
        

class DanilovEnvelopeMonitor22:
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
            "cov_00",
            "cov_01",
            "cov_02",
            "cov_03",
            "cov_11",
            "cov_12",
            "cov_13",
            "cov_22",
            "cov_23",
            "cov_33",
        ]:
            self.history[key] = []

    def package(self) -> None:
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
        
        params = env_params_from_bunch(bunch)
        param_matrix = env_vector_to_matrix(params)
        cov_matrix = 0.25 * np.matmul(param_matrix, param_matrix.T)

        self.history["s"].append(self.distance)
        self.history["xrms"].append(np.sqrt(cov_matrix[0, 0]))
        self.history["yrms"].append(np.sqrt(cov_matrix[2, 2]))

        for i in range(4):
            for j in range(i, 4):
                key = f"cov_{i}{j}"
                self.history[key] = cov_matrix[i, j]

        if self.verbose:
            message = ""
            message += "s={:0.3f} ".format(self.history["s"][-1])
            message += "xrms={:0.3f} ".format(self.history["xrms"][-1])
            message += "yrms={:0.3f} ".format(self.history["yrms"][-1])
            
    
class DanilovEnvelopeTracker22:
    def __init__(self, lattice: AccLattice, path_length_min: float= 1.00e-06) -> None:
        raise NotImplementedError