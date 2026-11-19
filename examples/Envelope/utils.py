import time
from typing import Any
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

from orbit.envelope import Envelope
from orbit.lattice import AccActionsContainer
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.utils.consts import mass_proton


def sample_distribution_gauss(size: int, dim: np.ndarray) -> np.ndarray:
    return np.random.multivariate_normal(
        mean=np.zeros(dim),
        cov=np.eye(dim),
        size=size,
    )


def sample_distribution_kv(size: int, dim: int) -> np.ndarray:
    X = np.random.normal(size=(size, dim))
    X /= np.linalg.norm(X, axis=1)[:, None]
    X /= np.std(X, axis=0)
    return X


def sample_distribution_waterbag(size: int, dim: int) -> np.ndarray:
    X = sample_distribution_kv(size, dim)
    dim = X.shape[1]
    r = np.random.uniform(size=size) ** (1.0 / dim)
    X *= r[:, None]
    X /= np.std(X, axis=0)
    return X


def sample_distribution(size: int, cov_matrix: np.ndarray, name: str) -> np.ndarray:
    dim = cov_matrix.shape[0]
    if name == "kv":
        X = sample_distribution_kv(size, dim)
    elif name == "waterbag":
        X = sample_distribution_waterbag(size, dim)
    elif name == "gauss":
        X = sample_distribution_gauss(size, dim)
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


def make_output_dir(script: str | Path, timestamp: bool = False) -> Path:
    output_dir = Path("outputs") / Path(script).stem
    if timestamp:
        output_dir = output_dir / time.strftime("%Y%m%d_%H%M%S")
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def make_fodo(
    kq: float = 0.25, nslice: int = 20, cells: int = 1
) -> tuple[TEAPOT_Lattice, list[QuadTEAPOT]]:
    lattice = TEAPOT_Lattice()
    quad_nodes = []
    for _ in range(cells):
        for length, strength in [
            (0.5, kq),
            (1.0, None),
            (1.0, -kq),
            (1.0, None),
            (0.5, kq),
        ]:
            node = (
                DriftTEAPOT(length=length)
                if strength is None
                else QuadTEAPOT(length=length, kq=strength)
            )
            node.setnParts(nslice)
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
            lattice.addNode(node)
            if strength is not None:
                quad_nodes.append(node)
    lattice.initialize()
    return lattice, quad_nodes


def make_envelope(
    parameters: np.ndarray | list[float],
    intensity: float = 5e9,
    emittance: float = 0.25e-6,
    kin_energy: float = 0.0025,
) -> Envelope:
    log_beta_x, alpha_x, log_beta_y, alpha_y = parameters
    beta_x = np.exp(log_beta_x)
    beta_y = np.exp(log_beta_y)

    cov_matrix = np.zeros((6, 6))
    cov_matrix[0, 0] = emittance * beta_x
    cov_matrix[2, 2] = emittance * beta_y
    cov_matrix[0, 1] = cov_matrix[1, 0] = -emittance * alpha_x
    cov_matrix[2, 3] = cov_matrix[3, 2] = -emittance * alpha_y
    cov_matrix[1, 1] = emittance * (1.0 + alpha_x**2) / beta_x
    cov_matrix[3, 3] = emittance * (1.0 + alpha_y**2) / beta_y
    cov_matrix[4, 4] = (5.0) ** 2 / 12.0
    return Envelope(
        cov_matrix=cov_matrix,
        kin_energy=kin_energy,
        mass=mass_proton,
        intensity=intensity,
    )


def twiss_parameters(envelope: Envelope) -> np.ndarray:
    result = []
    for i in (0, 2):
        c = envelope.cov_matrix[i : i + 2, i : i + 2]
        emittance = np.sqrt(np.linalg.det(c))
        beta = c[0, 0] / emittance
        alpha = -c[0, 1] / emittance
        result.extend([np.log(beta), alpha])
    return np.array(result)


def matched_envelope(
    lattice: AccLattice,
    intensity: float = 5e9,
    initial: np.ndarray | list[float] | None = None,
    optimizer_kwargs: dict = None,
) -> tuple[Envelope, float]:

    if initial is None:
        envelope = make_envelope([0, 0, 0, 0])
        matrix = lattice.trackEnvelope(envelope, fit=False, return_matrix=True)
        initial = []
        for i in (0, 2):
            m = matrix[i : i + 2, i : i + 2]
            cosine = np.trace(m) / 2
            if abs(cosine) >= 1:
                raise ValueError("Bare lattice is unstable")
            sine = np.copysign(np.sqrt(1 - cosine**2), m[0, 1])
            initial.extend([np.log(m[0, 1] / sine), (m[0, 0] - m[1, 1]) / (2 * sine)])

    def residual(parameters: np.ndarray) -> np.ndarray:
        envelope = make_envelope(parameters, intensity)
        lattice.trackEnvelope(envelope, fit=False, static=True, sc="2d")
        return twiss_parameters(envelope) - parameters

    optimizer_kwargs = optimizer_kwargs or {}
    optimizer_kwargs.setdefault("xtol", 1e-10)
    optimizer_kwargs.setdefault("ftol", 1e-10)
    optimizer_kwargs.setdefault("gtol", 1e-10)
    optimizer_kwargs.setdefault("verbose", 2)

    result = scipy.optimize.least_squares(residual, initial, **optimizer_kwargs)
    error = np.linalg.norm(result.fun)
    if not result.success or error > 1e-7:
        raise RuntimeError(f"Matching failed: residual={error:.3g}")
    return make_envelope(result.x, intensity), error


def track_history(
    lattice: AccLattice,
    envelope: Envelope,
    sc: str | None = None,
    static: bool = False,
) -> tuple[dict[str, np.ndarray], np.ndarray]:

    history = {}
    for key in ["s", "centroid", "covariance"]:
        history[key] = []

    def monitor(params_dict: dict) -> None:
        envelope = params_dict["envelope"]
        position = params_dict["path_length"]
        if history["s"] and history["s"][-1] == position:
            return
        history["s"].append(position)
        history["centroid"].append(envelope.centroid.copy())
        history["covariance"].append(envelope.cov_matrix.copy())

    action_container = AccActionsContainer()
    action_container.addAction(monitor, AccActionsContainer.ENTRANCE)
    action_container.addAction(monitor, AccActionsContainer.BODY)
    action_container.addAction(monitor, AccActionsContainer.EXIT)

    matrix = lattice.trackEnvelope(
        envelope,
        sc=sc,
        fit=False,
        static=static,
        action_container=action_container,
        return_matrix=True,
    )
    for key in history:
        history[key] = np.array(history[key])
    return history, matrix
