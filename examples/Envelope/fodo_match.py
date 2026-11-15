"""Compute periodic KV envelope in FODO lattice."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from utils import make_fodo
from utils import matched_envelope
from utils import make_output_dir
from utils import track_history
from utils import sample_distribution


def main(args: argparse.Namespace) -> None:
    output_dir = make_output_dir(__file__)
    np.random.seed(23)
    plt.style.use("style.mplstyle")

    lattice, _ = make_fodo(args.kq, args.nslice)
    envelope, residual = matched_envelope(lattice, args.intensity)

    initial = envelope.copy()
    history, _ = track_history(lattice, envelope.copy(), sc="2d", static=True)

    closure = []
    for _ in range(args.turns):
        lattice.trackEnvelope(envelope, sc="2d", fit=False, static=True)
        closure.append(
            np.linalg.norm(envelope.cov_matrix[:4, :4] - initial.cov_matrix[:4, :4])
            / np.linalg.norm(initial.cov_matrix[:4, :4])
        )
    print(f"Matching residual: {residual:.3g}; maximum relative closure error: {max(closure):.3g}")

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(history["s"], 1e3 * np.sqrt(history["covariance"][:, 0, 0]), label="x")
    ax.plot(history["s"], 1e3 * np.sqrt(history["covariance"][:, 2, 2]), label="y")
    ax.set(xlabel="s [m]", ylabel="RMS size [mm]")
    fig.savefig(output_dir / "matching.png")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kq", type=float, default=0.25)
    parser.add_argument("--intensity", type=float, default=5e9)
    parser.add_argument("--nslice", type=int, default=20)
    parser.add_argument("--nparts", type=int, default=10000)
    parser.add_argument("--turns", type=int, default=20)
    main(parser.parse_args())
