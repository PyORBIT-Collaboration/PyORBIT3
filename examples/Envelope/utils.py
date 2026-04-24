import numpy as np


def gen_dist_gauss(n: int, cov_matrix: np.ndarray) -> np.ndarray:
    return np.random.multivariate_normal(
        mean=np.zeros(cov_matrix.shape[0]),
        cov=cov_matrix,
        size=n,
    )


def gen_dist_kv(n: int, cov_matrix: np.ndarray) -> np.ndarray:
    X = np.random.normal(size=(n, cov_matrix.shape[0]))
    X /= np.linalg.norm(X, axis=1)[:, None]
    X /= np.std(X, axis=0)
    return X


def gen_dist_waterbag(n: int, cov_matrix: np.ndarray) -> np.ndarray:
    X = gen_dist_kv(n, cov_matrix)
    dim = X.shape[1]
    r = np.random.uniform(size=n) ** (1.0 / dim)
    X *= r[:, None]
    X /= np.std(X, axis=0)
    return X


def gen_dist(n: int, cov_matrix: np.ndarray, name: str) -> np.ndarray:
    if name == "kv":
        X = gen_dist_kv(n, cov_matrix)
    elif name == "waterbag":
        X = gen_dist_waterbag(n, cov_matrix)
    elif name == "gauss":
        X = gen_dist_gauss(n, cov_matrix)
    else:
        raise ValueError(f"Invalid distribution name: {name}")

    L = np.linalg.cholesky(cov_matrix)
    return np.matmul(X, L.T)


def build_rotation_matrix_xy(angle: float) -> np.ndarray:
    cs = np.cos(angle)
    sn = np.sin(angle)

    matrix = np.identity(4)
    matrix[0, 0] = matrix[1, 1] = +cs
    matrix[0, 2] = matrix[1, 3] = +sn
    matrix[2, 0] = matrix[3, 1] = -sn
    matrix[2, 2] = matrix[3, 3] = +cs
    return matrix


def project_cov_matrix(cov_matrix: np.ndarray, axis: tuple[int, ...]) -> np.ndarray:
    cov_matrix_proj = np.zeros((len(axis), len(axis)))
    for i in range(len(axis)):
        for j in range(len(axis)):
            cov_matrix_proj[i, j] = cov_matrix[axis[i], axis[j]]
    return cov_matrix_proj
