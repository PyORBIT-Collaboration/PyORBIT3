import numpy as np

from orbit.core.bunch import Bunch
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice


def build_rotation_matrix(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, s], [-s, c]])


def build_rotation_matrix_xy(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0, s, 0], [0, c, 0, s], [-s, 0, c, 0], [0, -s, 0, c]])


def build_phase_advance_matrix(*phase_advances: list[float]) -> np.ndarray:
    n = len(phase_advances)
    M = np.zeros((2 * n, 2 * n))
    for i, phase_advance in enumerate(phase_advances):
        i = i * 2
        M[i : i + 2, i : i + 2] = build_rotation_matrix(phase_advance)
    return M


def build_norm_matrix_from_twiss_2d(alpha: float, beta: float) -> np.ndarray:
    norm_matrix_inv = np.array([[beta, 0.0], [-alpha, 1.0]]) * np.sqrt(1.0 / beta)
    norm_matrix = np.linalg.inv(norm_matrix_inv)
    return norm_matrix


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
    classical_proton_radius = 1.53469e-18  # [m]
    gamma = 1.0 + (kin_energy / mass)  # Lorentz factor
    beta = np.sqrt(1.0 - (1.0 / gamma) ** 2)  # velocity/speed_of_light
    return (2.0 * classical_proton_radius * line_density) / (beta**2 * gamma**3)


def calc_cov_twiss(cov_matrix: np.ndarray) -> tuple[float, float, float]:
    emittance = np.sqrt(np.linalg.det(cov_matrix))
    alpha = -cov_matrix[0, 1] / emittance
    beta = cov_matrix[0, 0] / emittance
    return (alpha, beta, emittance)


def bunch_to_numpy(bunch: Bunch) -> np.ndarray:
    x = np.zeros((bunch.getSize(), 6))
    for i in range(bunch.getSize()):
        x[i, :] = [bunch.x(i), bunch.xp(i), bunch.y(i), bunch.yp(i), bunch.z(i), bunch.dE(i)]
    return x
