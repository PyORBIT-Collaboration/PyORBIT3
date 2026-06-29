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


def calc_norm_matrix_from_eigvecs(*eigvecs: list[np.ndarray]) -> np.ndarray:
    ndim = eigvecs[0].shape[0]

    V = np.zeros((ndim, ndim))
    for i, v in enumerate(eigvecs):
        V[:, 2 * i] = np.real(v)
        V[:, 2 * i + 1] = -np.imag(v)
    return np.linalg.inv(V)


def calc_norm_matrix_from_tmat(tmat: np.ndarray) -> np.ndarray:
    eig_res = np.linalg.eig(tmat)
    eigvecs = eig_res.eigenvectors[:, ::2].T
    for i in range(eigvecs.shape[0]):
        eigvecs[i] = normalize_eigvec(eigvecs[i])
    return calc_norm_matrix_from_eigvecs(eigvecs)


def calc_norm_matrix_from_cov(sigma: np.ndarray) -> np.ndarray:
    S = sigma
    U = build_poisson_matrix(S.shape[0])

    eig_res = np.linalg.eig(S @ U)
    eigvecs = eig_res.eigenvectors[:, ::2].T
    for i in range(eigvecs.shape[0]):
        eigvecs[i] = normalize_eigvec(eigvecs[i])
    return calc_norm_matrix_from_eigvecs(eigvecs)
