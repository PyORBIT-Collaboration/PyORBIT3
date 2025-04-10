import numpy as np


def unit_symplectic_matrix(ndim: int) -> np.ndarray:
    U = np.zeros((ndim, ndim))
    for i in range(0, ndim, 2):
        U[i: i + 2, i: i + 2] = [[0.0, 1.0], [-1.0, 0.0]]
    return U


def is_symplectic(M: np.ndarray, tol: float = 1.00e-06) -> bool:
    U = unit_symplectic_matrix(M.shape[0])
    residuals = U - np.linalg.multi_dot([M.T, U, M])
    return np.sum(residuals**2) < tol


def is_coupled(M: np.ndarray) -> bool:
    if M.shape[0] < 4:
        return False
    mask = np.zeros(M.shape)
    for i in range(0, M.shape[0], 2):
        mask[i: i + 2, i: i + 2] = 1.0
    return np.any(np.ma.masked_array(M, mask=mask))


def is_stable(M: np.ndarray, tol: float = 1.00e-08) -> np.ndarray:
    return all_eigvals_on_unit_circle(np.linalg.eigenvalues(M), tol=tol)


def all_eigvals_on_unit_circle(eigvals: np.ndarray, tol: float = 1.00e-08) -> bool:
    for eigval in eigvals:
        if abs(np.linalg.norm(eigval) - 1.0) > tol:
            return False
    return True


def eigtunes_from_eigvals(eigvals: np.ndarray) -> np.ndarray:
    def eigtune_from_eigval(eigval: np.complex) -> float:
        return np.arccos(eigval.real) / (2.0 * np.pi)
    
    return np.array([eigtune_from_eigval(eigvals[k]) for k in (0, 2)])


def phase_advance_matrix(*phase_advances) -> np.ndarray:
    def _rotation_matrix(angle: float) -> np.ndarray:
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[c, s], [-s, c]])

    ndim = 2 * len(phase_advances)
    M = np.zeros((ndim, ndim))
    for i in range(0, ndim, 2):
        M[i: i + 2, i: i + 2] = _rotation_matrix(phase_advances[i])
    return M


def cs_inv_norm_matrix_from_params(
    alpha_x: float,
    beta_x: float,
    alpha_y: float,
    beta_y: float,
) -> np.ndarray:
    def _build_matrix_2x2(alpha: float, beta: float) -> np.ndarray:
        return np.array([[beta, 0.0], [-alpha, 1.0]]) / np.sqrt(beta)

    V = np.zeros((4, 4))
    V[0:2, 0:2] = _build_matrix_2x2(alpha_x, beta_x)
    V[2:4, 2:4] = _build_matrix_2x2(alpha_y, beta_y)
    return V


def cs_params_from_transfer_matrix_2x2(M: np.ndarray) -> np.ndarray:
    params = dict()
    cos_phi = (M[0, 0] + M[1, 1]) / 2.0
    sign = 1.0
    if abs(M[0, 1]) != 0:
        sign = M[0, 1] / abs(M[0, 1])
    sin_phi = sign * np.sqrt(1.0 - cos_phi**2)
    beta = M[0, 1] / sin_phi
    alpha = (M[0, 0] - M[1, 1]) / (2.0 * sin_phi)
    params["alpha"] = alpha
    params["beta"] = beta
    # params["tune"] = np.arccos(cos_phi) / (2.0 * np.pi) * sign
    return params


def cs_params_from_transfer_matrix(M: np.ndarray) -> np.ndarray:
    params_x = cs_params_from_transfer_matrix_2x2(M[0:2, 0:2])
    params_y = cs_params_from_transfer_matrix_2x2(M[2:4, 2:4])
    return {
        "alpha_x": params_x["alpha"],
        "alpha_y": params_y["alpha"],
        "beta_x": params_x["beta"],
        "beta_y": params_y["beta"],
        # "tune_x": params_x["tune"],
        # "tune_y": params_y["tune"],
    }


def lb_normalize_eigvecs(eigvecs: np.ndarray) -> np.ndarray:
    ndim = eigvecs.shape[0]
    U = unit_symplectic_matrix(ndim)

    eigvecs = np.copy(eigvecs)
    for i in range(0, ndim, 2):
        v = eigvecs[:, i]
        val = np.linalg.multi_dot([np.conj(v), U, v]).imag
        if val > 0:
            eigvecs[:, i], eigvecs[:, i + 1] = eigvecs[:, i + 1], eigvecs[:, i]
        eigvecs[:, i : i + 2] *= np.sqrt(2 / np.abs(val))
    return eigvecs


def lb_inv_norm_matrix_from_eigvecs(eigvecs: np.ndarray) -> np.ndarray:
    V = np.zeros(eigvecs.shape)
    for i in range(0, V.shape[1], 2):
        V[:, i] = eigvecs[:, i].real
        V[:, i + 1] = (1j * eigvecs[:, i]).real
    return V


def lb_inv_norm_matrix_from_transfer_matrix(M: np.ndarray) -> np.ndarray:
    eigvals, eigvecs = np.linalg.eig(M)
    eigvecs = lb_normalize_eigvecs(eigvecs)
    return lb_inv_norm_matrix_from_eigvecs(eigvecs)


def lb_inv_norm_matrix_from_params_one_mode(
    alpha_lx: float, 
    beta_lx: float, 
    alpha_ly: float, 
    beta_ly: float, 
    u: float,
    nu: float,
    mode: int,
) -> np.ndarray:
    
    V = np.zeros((4, 4))
    if mode == 1:
        V[0, 0] = np.sqrt(beta_lx)
        V[0, 1] = 0.0
        V[1, 0] = -alpha_lx / np.sqrt(beta_lx)
        V[1, 1] = (1.0 - u) / np.sqrt(beta_lx)
        V[2, 0] = +np.sqrt(beta_ly) * np.cos(nu)
        V[2, 1] = -np.sqrt(beta_ly) * np.sin(nu)
        V[3, 0] = (u * np.sin(nu) - alpha_ly * np.cos(nu)) / np.sqrt(beta_ly)
        V[3, 1] = (u * np.cos(nu) + alpha_ly * np.sin(nu)) / np.sqrt(beta_ly)
    elif mode == 2:
        V[0, 2] = +np.sqrt(beta_lx) * np.cos(nu)
        V[0, 3] = -np.sqrt(beta_lx) * np.sin(nu)
        V[1, 2] = (u * np.sin(nu) - alpha_lx * np.cos(nu)) / np.sqrt(beta_lx)
        V[1, 3] = (u * np.cos(nu) + alpha_lx * np.sin(nu)) / np.sqrt(beta_lx)
        V[2, 2] = np.sqrt(beta_ly)
        V[2, 3] = 0.0
        V[3, 2] = -alpha_ly / np.sqrt(beta_ly)
        V[3, 3] = (1.0 - u) / np.sqrt(beta_ly)
    return V


def lb_params_from_transfer_matrix(M: np.ndarray) -> dict[str, float]:
    V = lb_inv_norm_matrix_from_eigvecs(M)

    beta_1x = V[0, 0] ** 2
    beta_2y = V[2, 2] ** 2
    alpha_1x = -np.sqrt(beta_1x) * V[1, 0]
    alpha_2y = -np.sqrt(beta_2y) * V[3, 2]
    u = 1.0 - (V[0, 0] * V[1, 1])
    nu1 = np.arctan2(-V[2, 1], V[2, 0])
    nu2 = np.arctan2(-V[0, 3], V[0, 2])
    beta_1y = (V[2, 0] / np.cos(nu1)) ** 2
    beta_2x = (V[0, 2] / np.cos(nu2)) ** 2
    alpha_1y = (u * np.sin(nu1) - V[3, 0] * np.sqrt(beta_1y)) / np.cos(nu1)
    alpha_2x = (u * np.sin(nu2) - V[1, 2] * np.sqrt(beta_2x)) / np.cos(nu2)
    return {
        "alpha_1x": alpha_1x,
        "alpha_1y": alpha_1y,
        "alpha_2x": alpha_2x,
        "alpha_2y": alpha_2y,
        "beta_1x": beta_1x,
        "beta_1y": beta_1y,
        "beta_2x": beta_2x,
        "beta_2y": beta_2y,
        "u": u,
        "nu1": nu1,
        "nu2": nu2,
    }