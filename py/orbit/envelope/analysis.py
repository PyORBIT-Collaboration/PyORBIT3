from typing import Any

import numpy as np

from .envelope import Envelope


class EnvelopeAnalysis:
    """Compute beam parameters from an envelope's current state.

    The envelope is retained by reference. Each method reads its current data,
    so an analysis object can be reused while the envelope is being tracked.
    Coordinates are ordered (x, x', y, y', z, dE).
    """

    def __init__(self, envelope: Envelope) -> None:
        self.envelope = envelope

    def centroid(self) -> np.ndarray:
        """Return a copy of the six centroid coordinates."""
        return self.envelope.centroid.copy()

    def covariance(self) -> np.ndarray:
        """Return a copy of the 6 x 6 covariance matrix."""
        return self.envelope.cov_matrix.copy()

    def rms(self, axis: int | None = None) -> np.ndarray | float:
        """Return rms sizes in x, y, z, or one plane (0, 1, 2)."""
        sizes = np.sqrt(np.diag(self.envelope.cov_matrix)[::2])
        return sizes if axis is None else sizes[axis]

    def emittances(self, normalized: bool = False) -> np.ndarray:
        """Return projected rms emittances for (x, x'), (y, y'), (z, dE).

        Normalization multiplies transverse values by beta * gamma and divides
        the longitudinal value by beta, following the envelope convention.
        """
        cov = self.envelope.cov_matrix
        values = np.array([
            np.sqrt(np.linalg.det(cov[i:i + 2, i:i + 2]))
            for i in (0, 2, 4)
        ])
        if normalized:
            beta = self.envelope.beta
            values *= np.array([beta * self.envelope.gamma] * 2 + [1.0 / beta])
        return values

    def eigenemittances(self) -> np.ndarray:
        """Return three sorted symplectic emittances in the envelope coordinates.

        Eigenvalues occur in conjugate pairs. Include zero modes for singular
        covariance matrices, such as a cold beam.
        """
        poisson = np.zeros((6, 6))
        for i in (0, 2, 4):
            poisson[i, i + 1] = 1.0
            poisson[i + 1, i] = -1.0
        eigenvalues = np.linalg.eigvals(self.envelope.cov_matrix @ poisson)
        return np.sort(np.abs(eigenvalues.imag))[::2]

    def to_dict(self) -> dict[str, Any]:
        """Return the standard beam parameters as an independent snapshot.

        For tracking, use an observer such as::

            def observe(envelope, s):
                return {"s": s, **EnvelopeAnalysis(envelope).to_dict()}
        """
        rms = self.rms()
        eps = self.emittances()
        eps_n = self.emittances(normalized=True)
        eigen = self.eigenemittances()
        return {
            "kin_energy": self.envelope.kin_energy,
            "gamma": self.envelope.gamma,
            "beta": self.envelope.beta,
            "mean": self.centroid(),
            "cov": self.covariance(),
            "rms_x": rms[0], "rms_y": rms[1], "rms_z": rms[2],
            "eps_x": eps[0], "eps_y": eps[1], "eps_z": eps[2],
            "eps_x_n": eps_n[0], "eps_y_n": eps_n[1], "eps_z_n": eps_n[2],
            "eps_1": eigen[0], "eps_2": eigen[1], "eps_3": eigen[2],
        }
