import math

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle


def get_dp_p_coeff(sync_part: SyncParticle) -> float:
    # dE/E = (beta^2) * dp/p
    # dE = (beta^2 * E) * dp/p
    # dE = (beta^2 * gamma * m * c^2) * dp/p
    beta = sync_part.beta()
    gamma = sync_part.gamma()
    rest_energy = sync_part.mass()  # GeV
    return 1.0 / (beta**2 * gamma * rest_energy)


def get_zp_coeff(sync_part: SyncParticle) -> float:
    # dE/E = (beta^2) * dp/p = (beta^2) * (gamma^2) z'
    # dE = (beta^2 * gamma^2 * E) * z'
    # dE = (beta^2 * gamma^3 * m * c^2) * z'
    beta = sync_part.beta()
    gamma = sync_part.gamma()
    rest_energy = sync_part.mass()
    return 1.0 / (beta**2 * gamma**3 * rest_energy)


def convert_matrix_dp_p_to_dE(matrix: np.ndarray, sync_part: SyncParticle) -> np.ndarray:
    # v = [x, x', y, y', z, dp/p]
    # w = [x, x', y, y', z, dE]
    # v = A w
    # v -> M v
    # w -> A M A^-1

    # scale = np.identity(7)
    # scale[5, 5] = dp_p_coeff
    #
    # scale_inv = np.identity(7)
    # scale_inv[5, 5] = 1.0 / dp_p_coeff
    #
    # return np.linalg.multi_dot([scale, matrix, scale_inv])

    dp_p_coeff = get_dp_p_coeff(sync_part)
    matrix[:5, 5] *= dp_p_coeff
    matrix[5, :5] /= dp_p_coeff
    matrix[5, 6] /= dp_p_coeff
    return matrix


def convert_matrix_zp_to_dE(matrix: np.ndarray, sync_part: SyncParticle) -> np.ndarray:
    zp_coeff = get_zp_coeff(sync_part)
    matrix[:5, 5] *= zp_coeff
    matrix[5, :5] /= zp_coeff
    matrix[5, 6] /= zp_coeff
    return matrix


def drift_matrix(length: float, sync_part: SyncParticle) -> np.ndarray:
    M = np.identity(7)
    M[0, 1] = length
    M[2, 3] = length
    M[4, 5] = length / sync_part.gamma() ** 2
    M = convert_matrix_dp_p_to_dE(M, sync_part)
    return M


def quad_matrix(length: float, kq: float, sync_part: SyncParticle) -> np.ndarray:
    if abs(kq) == 0:
        return drift_matrix(length=length, sync_part=sync_part)

    sqrt_abs_kq = math.sqrt(abs(kq))

    M = np.identity(7)
    if kq > 0:
        cx = np.cos(sqrt_abs_kq * length)
        sx = np.sin(sqrt_abs_kq * length)
        cy = np.cosh(sqrt_abs_kq * length)
        sy = np.sinh(sqrt_abs_kq * length)
        M[0, 0] = cx
        M[0, 1] = +sx / sqrt_abs_kq
        M[1, 0] = -sx * sqrt_abs_kq
        M[1, 1] = cx
        M[2, 2] = cy
        M[2, 3] = sy / sqrt_abs_kq
        M[3, 2] = sy * sqrt_abs_kq
        M[3, 3] = cy
    elif kq < 0:
        cx = np.cosh(sqrt_abs_kq * length)
        sx = np.sinh(sqrt_abs_kq * length)
        cy = np.cos(sqrt_abs_kq * length)
        sy = np.sin(sqrt_abs_kq * length)
        M[0, 0] = cx
        M[0, 1] = sx / sqrt_abs_kq
        M[1, 0] = sx * sqrt_abs_kq
        M[1, 1] = cx
        M[2, 2] = cy
        M[2, 3] = +sy / sqrt_abs_kq
        M[3, 2] = -sy * sqrt_abs_kq
        M[3, 3] = cy

    M[4, 5] = length / sync_part.gamma() ** 2
    M = convert_matrix_dp_p_to_dE(M, sync_part)
    return M


def bend_matrix(length: float, theta: float, sync_part: SyncParticle) -> np.ndarray:
    if length <= 0:
        return np.identity(7)

    rho = length / theta
    cx = math.cos(theta)
    sx = math.sin(theta)

    M = np.identity(7)
    M[0, 0] = cx
    M[0, 1] = rho * sx
    M[0, 5] = rho * (1.0 - cx)
    M[1, 0] = -sx / rho
    M[1, 1] = cx
    M[1, 5] = sx
    M[2, 3] = length
    M[4, 0] = -sx
    M[4, 1] = -rho * (1.0 - cx)
    M[4, 5] = -(sync_part.beta() ** 2) * length + rho * sx
    M = convert_matrix_dp_p_to_dE(M, sync_part)
    return M


def tilt_matrix(angle: float) -> np.ndarray:
    M = np.identity(7)
    M[0, 0] = M[1, 1] = +math.cos(angle)
    M[0, 2] = M[1, 3] = -math.sin(angle)
    M[2, 0] = M[3, 1] = +math.sin(angle)
    M[2, 2] = M[3, 3] = +math.cos(angle)
    return M


def translation_matrix(x: float = 0.0, y: float = 0.0, z: float = 0.0) -> np.ndarray:
    M = np.identity(7)
    M[0, -1] = x
    M[2, -1] = y
    M[4, -1] = z
    return M


def kick_matrix(kx: float = 0.0, ky: float = 0.0, kE: float = 0.0) -> np.ndarray:
    M = np.identity(7)
    M[1, -1] = kx
    M[3, -1] = ky
    M[5, -1] = kE
    return M


def solenoid_matrix(length: float, B: float, sync_part: SyncParticle) -> np.ndarray:
    if B == 0:
        return drift_matrix(length=length, sync_part=sync_part)

    phase = B * length

    V = np.identity(7)
    V[:4, :4] = 0.0
    V[0, 1] = -1.0 / B
    V[0, 2] = 0.5
    V[1, 0] = 0.5 * B
    V[1, 3] = 1.0
    V[2, 1] = 1.0 / B
    V[2, 2] = 0.5
    V[3, 0] = -0.5 * B
    V[3, 3] = 1.0

    M = np.identity(7)
    M[0, 0] = +1.0
    M[1, 1] = -1.0
    M[2, 2] = math.cos(phase)
    M[2, 3] = math.sin(phase) / B
    M[3, 2] = math.sin(phase) * (-B)
    M[3, 3] = math.cos(phase)
    M[4, 5] = length / sync_part.gamma() ** 2

    M = np.linalg.multi_dot([np.linalg.inv(V), M, V])
    M = convert_matrix_dp_p_to_dE(M, sync_part)
    return M


def cf_matrix(length: float, kq: float, sync_part: SyncParticle) -> np.ndarray:
    raise NotImplementedError()
