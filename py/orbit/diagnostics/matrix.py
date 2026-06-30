import numpy as np


def calc_eigtune(eigval: float) -> float:
    return np.arccos(np.real(eigval)) / (2.0 * np.pi)


def build_poisson_matrix(ndim: int) -> np.ndarray:
    U = np.zeros((ndim, ndim))
    for i in range(0, ndim, 2):
        U[i : i + 2, i : i + 2] = [[0.0, 1.0], [-1.0, 0.0]]
    return U


def normalize_eigvec(v: np.ndarray) -> np.ndarray:
    U = build_poisson_matrix(len(v))

    def complex_amplitude(v):
        return np.linalg.multi_dot([np.conj(v), U, v])

    if np.imag(complex_amplitude(v)) > 0.0:
        v = np.conj(v)

    v *= np.sqrt(2.0 / np.abs(complex_amplitude(v)))
    assert np.isclose(np.imag(complex_amplitude(v)), -2.0)
    assert np.isclose(np.real(complex_amplitude(v)), +0.0)
    return v


def build_norm_matrix_from_eigvecs(*eigvecs: list[np.ndarray]) -> np.ndarray:
    dim = len(eigvecs[0])

    V = np.zeros((dim, dim))
    for i, v in enumerate(eigvecs):
        V[:, 2 * i + 0] = +np.real(v)
        V[:, 2 * i + 1] = -np.imag(v)
    return np.linalg.inv(V)


def build_norm_matrix_from_tmat(matrix: np.ndarray) -> np.ndarray:
    eig_res = np.linalg.eig(matrix)
    eigvecs = eig_res.eigenvectors[:, ::2].T
    for i in range(eigvecs.shape[0]):
        eigvecs[i] = normalize_eigvec(eigvecs[i])
    return build_norm_matrix_from_eigvecs(*eigvecs)


def build_norm_matrix_from_cov(cov_matrix: np.ndarray) -> np.ndarray:
    S = cov_matrix
    U = build_poisson_matrix(S.shape[0])

    eig_res = np.linalg.eig(S @ U)
    eigvecs = eig_res.eigenvectors[:, ::2].T
    for i in range(eigvecs.shape[0]):
        eigvecs[i] = normalize_eigvec(eigvecs[i])
    return build_norm_matrix_from_eigvecs(*eigvecs)


class TransferMatrixAnalysis:
    def __init__(self, M: np.ndarray) -> None:
        self.M = M
        self.ndim = M.shape[0]

        self.eigvals, self.eigvecs = np.linalg.eig(M)
        self.eigvecs = self.eigvecs.T
        self.eigvecs = self.eigvecs[::2]
        for i in range(self.eigvecs.shape[0]):
            self.eigvecs[i] = normalize_eigvec(self.eigvecs[i])

        self.eigvals = self.eigvals[::2]
        self.eigtunes = [calc_eigtune(eigval) for eigval in self.eigvals]

        self.V_inv = build_norm_matrix_from_eigvecs(*self.eigvecs)
        self.V = np.linalg.inv(self.V_inv)

    def cov_matrix(self, *emittances: float) -> np.ndarray:
        S = np.eye(self.ndim)
        if len(emittances) > 0:
            S = np.diag(np.repeat(emittances, 2))
        return self.V @ S @ self.V.T
