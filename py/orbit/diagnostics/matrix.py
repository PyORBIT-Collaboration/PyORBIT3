import numpy as np


def calc_eigtune(eigval: float) -> float:
    return np.arccos(np.real(eigval)) / (2.0 * np.pi)


def build_poisson_matrix(ndim: int) -> np.ndarray:
    U = np.zeros((ndim, ndim))
    for i in range(0, ndim, 2):
        U[i : i + 2, i : i + 2] = [[0.0, 1.0], [-1.0, 0.0]]
    return U


def normalize_eigvec(v: np.ndarray) -> np.ndarray:
    v = np.copy(v)
    U = build_poisson_matrix(len(v))

    def complex_amplitude(v):
        return np.conj(v).T @ U @ v

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


def symplectic_eig(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eigvals, eigvecs = np.linalg.eig(matrix)
    eigvecs = eigvecs.T

    eigvals = eigvals[::2]
    eigvecs = eigvecs[::2]

    for i in range(eigvecs.shape[0]):
        eigvecs[i] = normalize_eigvec(eigvecs[i])

    return eigvals, eigvecs


def build_norm_matrix_from_tmat(matrix: np.ndarray) -> np.ndarray:
    eigvals, eigvecs = symplectic_eig(matrix)
    return build_norm_matrix_from_eigvecs(*eigvecs)


def build_norm_matrix_from_cov_uncoupled(cov_matrix: np.ndarray) -> np.ndarray:
    S = cov_matrix
    U = build_poisson_matrix(S.shape[0])

    V_inv = np.zeros_like(S)
    for i in range(0, S.shape[0], 2):
        eigvals, eigvecs = symplectic_eig(S[i:i+2, i:i+2] @ U[i:i+2, i:i+2])
        V_inv[i:i+2, i:i+2] = build_norm_matrix_from_eigvecs(*eigvecs)
    return V_inv


def build_norm_matrix_from_cov(cov_matrix: np.ndarray) -> np.ndarray:
    if not is_coupled(cov_matrix):
        return build_norm_matrix_from_cov_uncoupled(cov_matrix)

    S = cov_matrix
    U = build_poisson_matrix(S.shape[0])

    eigvals, eigvecs = symplectic_eig(S @ U)

    emittances = np.abs(np.imag(eigvals))
    if np.all(np.abs(emittances - emittances[0]) < 1e-15):
        raise ValueError("Eigenemittances are equal eigenvectors are degenerate and V will not be correct.")

    # Order by emittances
    order = np.argsort(emittances)
    eigvals = eigvals[order]
    eigvecs = eigvecs[order]
    return build_norm_matrix_from_eigvecs(*eigvecs)


def is_coupled(matrix: np.ndarray) -> bool:
    assert matrix.ndim == 2
    assert matrix.shape[0] == matrix.shape[1]
    assert matrix.shape[0] % 2 == 0
    matrix = np.copy(matrix)
    for i in range(0, matrix.shape[0], 2):
        matrix[i:i+2, i:i+2] = 0.0
    return not np.all(np.isclose(matrix, 0.0))


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

        self.is_coupled = is_coupled(M)

    def cov_matrix(self, *emittances: float) -> np.ndarray:
        S = np.eye(self.ndim)
        if len(emittances) > 0:
            S = np.diag(np.repeat(emittances, 2))
        return self.V @ S @ self.V.T
