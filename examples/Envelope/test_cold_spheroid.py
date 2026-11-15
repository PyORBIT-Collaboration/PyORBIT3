"""Expansion of cold uniform spheroid in free space."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.special

from orbit.envelope import Envelope
from orbit.teapot import DriftTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.utils.consts import mass_proton

from utils import make_output_dir
from utils import track_history


def main(args: argparse.Namespace) -> None:
    output_dir = make_output_dir(__file__)
    plt.style.use("style.mplstyle")

    radii = np.array([args.radius, args.radius, args.radius * args.aspect])
    envelope = Envelope(
        kin_energy=args.kin_energy,
        mass=mass_proton,
        intensity=args.intensity,
        cov_matrix=np.zeros((6, 6)),
    )

    # For a uniform ellipsoid, each spatial variance is semiaxis**2 / 5.
    # Longitudinal covariance is stored in lab coordinates, z_lab = z_rest / gamma.
    envelope.cov_matrix[0, 0] = radii[0] ** 2 / 5.0
    envelope.cov_matrix[2, 2] = radii[1] ** 2 / 5.0
    envelope.cov_matrix[4, 4] = (radii[2] / envelope.gamma) ** 2 / 5.0
    coefficient = (
        args.intensity
        * envelope.classical_radius
        / (envelope.beta * envelope.gamma) ** 2
    )

    def derivative(s: float, state: np.ndarray) -> np.ndarray:
        a, b, c = state[:3] ** 2
        acceleration = (
            coefficient
            * state[:3]
            * [
                scipy.special.elliprd(b, c, a),
                scipy.special.elliprd(a, c, b),
                scipy.special.elliprd(a, b, c),
            ]
        )
        return np.r_[state[3:], acceleration]

    lattice = TEAPOT_Lattice()
    for _ in range(args.nslice):
        lattice.addNode(DriftTEAPOT(length=args.length / args.nslice))
    lattice.initialize()

    gamma = envelope.gamma
    history, _ = track_history(lattice, envelope, sc="3d")

    tracked = np.zeros((history["covariance"].shape[0], 3))
    tracked[:, 0] = np.sqrt(5.0 * history["covariance"][:, 0, 0])
    tracked[:, 1] = np.sqrt(5.0 * history["covariance"][:, 2, 2])
    tracked[:, 2] = np.sqrt(5.0 * history["covariance"][:, 4, 4]) * gamma

    reference = scipy.integrate.solve_ivp(
        derivative,
        [0.0, args.length],
        np.r_[radii, np.zeros(3)],
        dense_output=True,
        rtol=1e-10,
        atol=1e-12,
    )
    if not reference.success:
        raise RuntimeError(reference.message)

    exact = reference.sol(history["s"])[:3].T
    print(f"Maximum relative radius error: {np.max(abs(tracked / exact - 1)):.3g}")

    fig, axs = plt.subplots(ncols=3, sharex=True, sharey=False, figsize=(10.0, 2.5))
    for i, dim in enumerate(("x", "y", "z")):
        axs[i].plot(history["s"], 1e3 * tracked[:, i], color="black", label=f"ENV")
        axs[i].plot(history["s"], 1e3 * exact[:, i], ls="--", color="red", label=f"ODE")
        axs[i].set_xlabel("s [m]")
        axs[i].set_ylabel(f"{dim} [mm]")
        axs[i].legend(loc="upper left", fontsize="small")
    fig.savefig(output_dir / "expansion")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=float, default=0.01)
    parser.add_argument("--aspect", type=float, default=2.0)
    parser.add_argument("--length", type=float, default=10.0)
    parser.add_argument("--kin-energy", type=float, default=0.01)
    parser.add_argument("--intensity", type=float, default=1e8)
    parser.add_argument("--nslice", type=int, default=400)
    main(parser.parse_args())
