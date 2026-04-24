import numpy as np

def gen_kv(n: int, cov_matrix: np.ndarray) -> np.ndarray:    
    rng = np.random.default_rng()
    X = rng.normal(size=(n, cov_matrix.shape[0]))
    X /= np.linalg.norm(X, axis=1)[:, None]
    X /= np.std(X, axis=0)
    X = np.matmul(X, np.linalg.cholesky(cov_matrix).T)
    return X