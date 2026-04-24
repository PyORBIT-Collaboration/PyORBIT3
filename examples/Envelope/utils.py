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