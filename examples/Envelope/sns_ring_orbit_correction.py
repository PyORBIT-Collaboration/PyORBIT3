"""Correct horizontal closed orbit based on centroid tracking."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

from orbit.envelope import Envelope
from orbit.lattice import AccNode
from orbit.teapot import KickTEAPOT
from orbit.teapot import TEAPOT_Ring
from orbit.utils.consts import mass_proton

from utils import make_output_dir
from utils import track_history


def main(args: argparse.Namespace) -> None:
    output_dir = make_output_dir(__file__)
    plt.style.use("style.mplstyle")

    lattice = TEAPOT_Ring()
    lattice.readMADX(str(Path(__file__).with_name("inputs") / "sns_ring.lat"), "rnginjsol")
    for node in lattice.getNodes():
        if hasattr(node, "setUsageFringeFieldIN"):
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
        if node.getLength() > 1.0:
            node.setnParts(1 + int(node.getLength()))

    nodes = lattice.getNodes()
    if not 1 <= args.correctors <= len(nodes) // 2:
        raise ValueError("Corrector count must be between 1 and half the ring node count")

    positions = np.linspace(0, len(nodes), 2 * args.correctors, endpoint=False, dtype=int)
    rng = np.random.default_rng(23)
    corrector_nodes = []
    for i, index in enumerate(positions):
        corrector_node = KickTEAPOT(
            name=f"orbit_kick_{i}", kx=(rng.normal(0.0, args.error) if i % 2 else 0.0)
        )
        nodes[index].addChildNode(corrector_node, AccNode.ENTRANCE)
        if i % 2 == 0:
            corrector_nodes.append(corrector_node)
    lattice.initialize()

    def closed_orbit() -> tuple[np.ndarray, np.ndarray, float]:
        beam = Envelope(kin_energy=1.3, mass=mass_proton)
        matrix = lattice.trackEnvelope(beam.copy(), fit=False, return_matrix=True)
        beam.centroid[:4] = np.linalg.solve(np.eye(4) - matrix[:4, :4], matrix[:4, 6])
        initial = beam.centroid.copy()
        history, _ = track_history(lattice, beam)
        s = history["s"]
        centroid = history["centroid"]
        closure = np.linalg.norm(centroid[-1, :4] - initial[:4])
        return s, centroid[:, 0], closure

    s, before, _ = closed_orbit()

    columns = []
    for corrector_node in corrector_nodes:
        corrector_node.setParam("kx", args.delta)
        _, displaced, _ = closed_orbit()
        columns.append((displaced - before) / args.delta)
        corrector_node.setParam("kx", 0.0)

    response = np.column_stack(columns)
    result = scipy.optimize.lsq_linear(
        response, -before, bounds=(-args.max_kick, args.max_kick), tol=1e-12
    )
    if not result.success:
        raise RuntimeError(result.message)

    for corrector_node, strength in zip(corrector_nodes, result.x):
        corrector_node.setParam("kx", strength)

    _, after, closure = closed_orbit()

    print(f"Orbit RMS: {1e3 * np.std(before):.4g} -> {1e3 * np.std(after):.4g} mm")
    print(
        f"Orbit norm/sqrt(BPMs): {1e3 * np.sqrt(np.mean(before**2)):.4g} -> "
        f"{1e3 * np.sqrt(np.mean(after**2)):.4g} mm; transverse closure: {closure:.3g}"
    )
    print("Corrector kicks [rad]:", result.x)

    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(s, 1e3 * before, color="black", ls=":", label="Before")
    ax.plot(s, 1e3 * after, color="black", ls="-", label="After")
    ax.set(xlabel="s [m]", ylabel="Horizontal closed orbit [mm]")
    ax.legend()
    fig.savefig(output_dir / "orbit.png")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--correctors", type=int, default=12)
    parser.add_argument("--error", type=float, default=1e-4, help="RMS injected dipole error [rad]")
    parser.add_argument("--max-kick", type=float, default=1e-3)
    parser.add_argument("--delta", type=float, default=1e-6)
    main(parser.parse_args())
