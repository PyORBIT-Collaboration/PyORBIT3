import numpy as np


def norm_matrix_from_twiss_cs(
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
    return np.linalg.inv(V)


def norm_matrix_from_twiss_lb_one_mode(
    alpha_lx: float, 
    alpha_ly: float, 
    beta_lx: float, 
    beta_ly: float, 
    nu: float,
    u: float,
    mode: int,
) -> np.ndarray:
    
    _cos = np.cos(nu)
    _sin = np.sin(nu)
    
    V = np.zeros((4, 4))
    if mode == 1:
        V[0, 0] = np.sqrt(beta_lx)
        V[0, 1] = 0.0
        V[1, 0] = -alpha_lx / np.sqrt(beta_lx)
        V[1, 1] = (1.0 - u) / np.sqrt(beta_lx)
        V[2, 0] = +np.sqrt(beta_ly) * _cos
        V[2, 1] = -np.sqrt(beta_ly) * _sin
        V[3, 0] = (u * _sin - alpha_ly * _cos) / np.sqrt(beta_ly)
        V[3, 1] = (u * _cos + alpha_ly * _sin) / np.sqrt(beta_ly)
    else:
        V[0, 2] = +np.sqrt(beta_lx) * _cos
        V[0, 3] = -np.sqrt(beta_lx) * _sin
        V[1, 2] = (u * _sin - alpha_lx * _cos) / np.sqrt(beta_lx)
        V[1, 3] = (u * _cos + alpha_lx * _sin) / np.sqrt(beta_lx)
        V[2, 2] = np.sqrt(beta_ly)
        V[2, 3] = 0.0
        V[3, 2] = -alpha_ly / np.sqrt(beta_ly)
        V[3, 3] = (1.0 - u) / np.sqrt(beta_ly)
    return np.linalg.inv(V)
