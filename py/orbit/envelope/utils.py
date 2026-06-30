import numpy as np

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
    dp_p_coeff = get_dp_p_coeff(sync_part)
    matrix[:5, 5] *= dp_p_coeff
    matrix[5, :5] /= dp_p_coeff
    matrix[5, 6] /= dp_p_coeff  # driving term
    return matrix


def convert_matrix_zp_to_dE(matrix: np.ndarray, sync_part: SyncParticle) -> np.ndarray:
    zp_coeff = get_zp_coeff(sync_part)
    matrix[:5, 5] *= zp_coeff
    matrix[5, :5] /= zp_coeff
    matrix[5, 6] /= zp_coeff  # driving term
    return matrix


def gen_dist_gauss(size: int, cov_matrix: np.ndarray) -> np.ndarray:
    return np.random.multivariate_normal(
        mean=np.zeros(cov_matrix.shape[0]),
        cov=cov_matrix,
        size=size,
    )


def gen_dist_kv(size: int, cov_matrix: np.ndarray) -> np.ndarray:
    X = np.random.normal(size=(size, cov_matrix.shape[0]))
    X /= np.linalg.norm(X, axis=1)[:, None]
    X /= np.std(X, axis=0)
    return X


def gen_dist_waterbag(size: int, cov_matrix: np.ndarray) -> np.ndarray:
    X = gen_dist_kv(size, cov_matrix)
    dim = X.shape[1]
    r = np.random.uniform(size=size) ** (1.0 / dim)
    X *= r[:, None]
    X /= np.std(X, axis=0)
    return X


def gen_dist(size: int, cov_matrix: np.ndarray, name: str) -> np.ndarray:
    if name == "kv":
        X = gen_dist_kv(size, cov_matrix)
    elif name == "waterbag":
        X = gen_dist_waterbag(size, cov_matrix)
    elif name == "gauss":
        X = gen_dist_gauss(size, cov_matrix)
    else:
        raise ValueError(f"Invalid distribution name: {name}")

    L = np.linalg.cholesky(cov_matrix)
    return np.matmul(X, L.T)


def proj_cov_matrix(cov_matrix: np.ndarray, axis: tuple[int, ...]) -> np.ndarray:
    cov_matrix_proj = np.zeros((len(axis), len(axis)))
    for i in range(len(axis)):
        for j in range(len(axis)):
            cov_matrix_proj[i, j] = cov_matrix[axis[i], axis[j]]
    return cov_matrix_proj
