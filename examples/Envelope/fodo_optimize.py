"""Optimize quadrupoles to minimize RMS size fluctuations in FODO lattice."""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

from utils import make_fodo
from utils import make_output_dir
from utils import matched_envelope
from utils import track_history
from utils import twiss_parameters


def main(args: argparse.Namespace) -> None:
    output_dir = make_output_dir(__file__)
    plt.style.use("style.mplstyle")

    cell, _ = make_fodo(args.kq, args.nslice)
    envelope, _ = matched_envelope(cell, args.intensity)
    lattice, quads = make_fodo(args.kq, args.nslice, args.cells)
    signs = np.sign([node.getParam("kq") for node in quads])
    initial = args.kq * np.random.default_rng(23).uniform(0.7, 1.3, len(quads))
    size_scale = np.sqrt(envelope.cov_matrix[[0, 2], [0, 2]])
    target = twiss_parameters(envelope)

    def evaluate(strengths: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        for node, strength, sign in zip(quads, strengths, signs):
            node.setParam("kq", sign * strength)
        envelope_out = envelope.copy()
        history, _ = track_history(lattice, envelope_out, sc="2d")
        sizes = np.sqrt(history["covariance"][:, [0, 2], [0, 2]])
        return history["s"], sizes, twiss_parameters(envelope_out)

    def residual(strengths: np.ndarray) -> np.ndarray:
        _, sizes, endpoint = evaluate(strengths)
        fluctuations = (sizes - sizes.mean(axis=0)) / size_scale
        return np.r_[
            fluctuations.ravel() / np.sqrt(len(sizes)),
            args.match_weight * (endpoint - target),
        ]

    def objective(strengths: np.ndarray) -> float:
        return 0.5 * np.mean(np.abs(residual(strengths)))

    s, before, endpoint_before = evaluate(initial)
    result = scipy.optimize.least_squares(
        residual,
        initial,
        bounds=(0.05 * args.kq, 3 * args.kq),
        max_nfev=args.max_nfev,
        ftol=1e-12,
        gtol=1e-12,
        verbose=2,
    )
    _, after, endpoint_after = evaluate(result.x)

    print(result.message)
    for label, sizes, endpoint in [
        ("Before", before, endpoint_before),
        ("After", after, endpoint_after),
    ]:
        print(
            f"{label}: normalized size fluctuation={np.linalg.norm(np.std(sizes / size_scale, axis=0)):.4g}; "
            f"endpoint Twiss residual={np.linalg.norm(endpoint - target):.4g}"
        )
    print("Signed quadrupole strengths [1/m²]:", signs * result.x)

    fig, ax = plt.subplots()
    for i, name in enumerate(("x", "y")):
        (line,) = ax.plot(s, 1e3 * after[:, i], label=f"{name}: optimized")
        ax.plot(
            s,
            1e3 * before[:, i],
            "--",
            color=line.get_color(),
            label=f"{name}: initial",
        )
    ax.set(xlabel="s [m]", ylabel="RMS size [mm]")
    ax.legend()
    fig.savefig(output_dir / "optimization.png")
    plt.close(fig)
    print(output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=4)
    parser.add_argument("--kq", type=float, default=0.25)
    parser.add_argument("--intensity", type=float, default=5e9)
    parser.add_argument("--nslice", type=int, default=10)
    parser.add_argument("--match-weight", type=float, default=1.0)
    parser.add_argument("--max-nfev", type=int, default=120)
    main(parser.parse_args())
