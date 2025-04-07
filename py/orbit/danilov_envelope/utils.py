import numpy as np

from orbit.core.bunch import Bunch
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice


def rotation_matrix(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, s], [-s, c]])


def rotation_matrix_xy(angle: float) -> np.ndarray:
    """4 x 4 matrix to rotate [x, x', y, y'] clockwise in the x-y plane (angle in radians)."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0, s, 0], [0, c, 0, s], [-s, 0, c, 0], [0, -s, 0, c]])


def phase_advance_matrix(*phase_advances) -> np.ndarray:
    """Phase advance matrix (clockwise rotation in each phase plane).

    Parameters
    ---------
    mu1, mu2, ..., mun : float
        The phase advance in each plane.

    Returns
    -------
    ndarray, shape (2n, 2n)
        Matrix which rotates x-x', y-y', z-z', etc. by the phase advances.
    """
    n = len(phase_advances)
    M = np.zeros((2 * n, 2 * n))
    for i, phase_advance in enumerate(phase_advances):
        i = i * 2
        M[i : i + 2, i : i + 2] = rotation_matrix(phase_advance)
    return M


def get_transfer_matrix(lattice: AccLattice, bunch: Bunch) -> np.ndarray:
    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
    M = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            M[i, j] = matrix_lattice.oneTurnMatrix.get(i, j)
    return M


def fit_transfer_matrix(lattice: AccLattice, bunch: Bunch) -> np.ndarray:
    step_arr_init = np.full(4, 1.00e-06)
    step_arr = np.copy(step_arr_init)
    step_reduce = 20.0

    _bunch = Bunch()
    bunch.copyEmptyBunchTo(_bunch)

    _bunch.addParticle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    _bunch.addParticle(step_arr[0] / step_reduce, 0.0, 0.0, 0.0, 0.0, 0.0)
    _bunch.addParticle(0.0, step_arr[1] / step_reduce, 0.0, 0.0, 0.0, 0.0)
    _bunch.addParticle(0.0, 0.0, step_arr[2] / step_reduce, 0.0, 0.0, 0.0)
    _bunch.addParticle(0.0, 0.0, 0.0, step_arr[3] / step_reduce, 0.0, 0.0)
    _bunch.addParticle(step_arr[0], 0.0, 0.0, 0.0, 0.0, 0.0)
    _bunch.addParticle(0.0, step_arr[1], 0.0, 0.0, 0.0, 0.0)
    _bunch.addParticle(0.0, 0.0, step_arr[2], 0.0, 0.0, 0.0)
    _bunch.addParticle(0.0, 0.0, 0.0, step_arr[3], 0.0, 0.0)

    lattice.trackBunch(_bunch)

    X = get_bunch_coords(bunch)
    X = X[:, (0, 2, 3, 4)]

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


def get_perveance(mass: float, kin_energy: float, line_density: float) -> float:
    """Compute dimensionless space charge perveance.

    Parameters
    ----------
    mass : float
        Mass per particle [GeV/c^2].
    kin_energy : float
        Kinetic energy per particle [GeV].
    line_density : float
        Number density in longitudinal plane [m^-1].

    Returns
    -------
    float
        Dimensionless space charge perveance.
    """
    classical_proton_radius = 1.53469e-18  # [m]
    gamma = 1.0 + (kin_energy / mass)  # Lorentz factor
    beta = np.sqrt(1.0 - (1.0 / gamma) ** 2)  # velocity/speed_of_light
    return (2.0 * classical_proton_radius * line_density) / (beta**2 * gamma**3)


def calc_twiss_2d(cov_matrix: np.ndarray) -> tuple[float, float, float]:
    emittance = np.sqrt(np.linalg.det(cov_matrix))
    alpha = -cov_matrix[0, 1] / emittance
    beta = cov_matrix[0, 0] / emittance
    return (alpha, beta, emittance)


def get_bunch_coords(bunch: Bunch) -> np.ndarray:
    x = np.zeros((bunch.getSize(), 6))
    for i in range(bunch.getSize()):
        x[i, :] = [bunch.x(i), bunch.xp(i), bunch.y(i), bunch.yp(i), bunch.z(i), bunch.dE(i)]
    return x