"""Test envelope tracker in SNS ring."""

import argparse
import copy
import math
import os
import pathlib
import time

import numpy as np
import matplotlib.pyplot as plt

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.bunch_utils import collect_bunch
from orbit.envelope import Envelope
from orbit.envelope import EnvelopeTracker
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_Ring
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from plot import plot_rms_ellipse
from plot import plot_corner
from utils import gen_dist
from utils import build_rotation_matrix_xy
from utils import project_cov_matrix

plt.style.use("style.mplstyle")


def main(args: argparse.Namespace) -> None:

    # Setup
    # ------------------------------------------------------------------------------

    path = pathlib.Path(__file__)
    output_dir = os.path.join("outputs", path.stem)
    os.makedirs(output_dir, exist_ok=True)

    # Create lattice
    # ------------------------------------------------------------------------------

    lattice = TEAPOT_Lattice()
    lattice.readMADX("inputs/sns_ring.lat", "rnginjsol")
    lattice.initialize()

    for node in lattice.getNodes():
        try:
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
        except:
            pass

    if args.sol:
        for name in ["scbdsol_c13a", "scbdsol_c13b"]:
            node = lattice.getNodeForName(name)
            node.setParam("B", 0.15)

    for node in lattice.getNodes():
        max_length = 1.0
        if node.getLength() > max_length:
            node.setnParts(1 + int(node.getLength() / max_length))

    # Create envelope
    # ------------------------------------------------------------------------------

    # Create bunch
    bunch = Bunch()
    bunch.mass(mass_proton)
    sync_part = bunch.getSyncParticle()
    sync_part.kinEnergy(args.kin_energy)

    # Find periodic lattice parameters
    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
    matrix_lattice_params = matrix_lattice.getRingParametersDict()
    alpha_x = matrix_lattice_params["alpha x"]
    alpha_y = matrix_lattice_params["alpha y"]
    beta_x = matrix_lattice_params["beta x [m]"]
    beta_y = matrix_lattice_params["beta y [m]"]
    eps_x = 25.0e-06
    eps_y = eps_x

    print(matrix_lattice_params)

    # Generate covariance matrix
    cov_matrix = np.zeros((6, 6))
    cov_matrix[0, 0] = eps_x * beta_x
    cov_matrix[2, 2] = eps_y * beta_y
    cov_matrix[0, 1] = cov_matrix[1, 0] = -eps_x * alpha_x
    cov_matrix[2, 3] = cov_matrix[3, 2] = -eps_y * alpha_y
    cov_matrix[1, 1] = eps_x * (1.0 + alpha_x**2) / beta_x
    cov_matrix[3, 3] = eps_y * (1.0 + alpha_y**2) / beta_y
    cov_matrix[4, 4] = (args.bunch_length / 4.0) ** 2
    cov_matrix[5, 5] = 0.0

    # Tilt
    if args.tilt:
        rot_matrix = np.identity(6)
        rot_matrix[:4, :4] = build_rotation_matrix_xy(angle=(args.tilt * math.pi))
        cov_matrix = np.linalg.multi_dot([rot_matrix, cov_matrix, rot_matrix.T])

    # Mismatch
    cov_matrix[0, 0] *= (1.0 + args.mismatch_x) ** 2
    cov_matrix[2, 2] *= (1.0 + args.mismatch_y) ** 2
    cov_matrix_init = np.copy(cov_matrix)

    # Offset
    centroid_init = np.zeros(6)
    centroid_init[0] += args.offset_x
    centroid_init[2] += args.offset_y

    # Create envelope
    envelope = Envelope(
        bunch=bunch,
        cov_matrix=cov_matrix_init,
        centroid=centroid_init,
        intensity=args.intensity,
    )

    # Track envelope
    # ------------------------------------------------------------------------------

    print("TRACK ENVELOPE")

    tracker = EnvelopeTracker(
        lattice,
        handle_unknown=args.handle_unknown,
        space_charge=("2d" if args.sc else None),
    )

    history = {"xrms": [], "yrms": [], "xavg": [], "yavg": []}
    start_time = time.time()

    for turn in range(args.turns):
        if turn > 0:
            tracker.track(envelope)

        cov_matrix = envelope.cov()
        centroid = envelope.centroid()

        xrms = 1000.0 * math.sqrt(cov_matrix[0, 0])
        yrms = 1000.0 * math.sqrt(cov_matrix[2, 2])
        xavg = 1000.0 * centroid[0]
        yavg = 1000.0 * centroid[2]

        message = ""
        message += " turn={}".format(turn)
        message += " time={:0.2f}".format(time.time() - start_time)
        message += " xrms={:0.2f}".format(xrms)
        message += " yrms={:0.2f}".format(yrms)
        message += " xavg={:0.2f}".format(xavg)
        message += " yavg={:0.2f}".format(yavg)
        print(message)

        history["xrms"].append(xrms)
        history["yrms"].append(yrms)
        history["xavg"].append(xavg)
        history["yavg"].append(yavg)

    histories = {}
    histories["envelope"] = copy.deepcopy(history)

    # Track bunch
    # ------------------------------------------------------------------------------

    print("TRACK BUNCH")

    rng = np.random.default_rng()

    bunch_coords = np.zeros((args.nparts, 6))
    bunch_coords[:, :4] = gen_dist(
        size=args.nparts, cov_matrix=cov_matrix_init[0:4, 0:4], name=args.dist
    )
    bunch_coords[:, 4] = args.bunch_length * rng.uniform(-0.5, 0.5, size=args.nparts)
    bunch_coords += centroid_init[None, :6]

    for i in range(bunch_coords.shape[0]):
        bunch.addParticle(*bunch_coords[i])

    if args.sc:
        sc_calc = SpaceChargeCalc2p5D(64, 64, 1)
        sc_path_length_min = 1.00e-06
        sc_nodes = setSC2p5DAccNodes(lattice, sc_path_length_min, sc_calc)

        bunch_size = bunch.getSizeGlobal()
        bunch.macroSize(args.intensity / bunch_size)

    history = {"xrms": [], "yrms": [], "xavg": [], "yavg": []}
    start_time = time.time()

    for turn in range(args.turns):
        if turn > 0:
            lattice.trackBunch(bunch)

        twiss_calc = BunchTwissAnalysis()
        twiss_calc.computeBunchMoments(bunch, 2, 0, 0)

        cov_matrix = np.zeros((6, 6))
        for i in range(6):
            for j in range(i + 1):
                cov_matrix[i, j] = twiss_calc.getCorrelation(j, i)
                cov_matrix[j, i] = cov_matrix[i, j]

        xrms = 1000.0 * math.sqrt(cov_matrix[0, 0])
        yrms = 1000.0 * math.sqrt(cov_matrix[2, 2])
        xavg = 1000.0 * twiss_calc.getAverage(0)
        yavg = 1000.0 * twiss_calc.getAverage(2)

        message = ""
        message += " turn={}".format(turn)
        message += " time={:0.2f}".format(time.time() - start_time)
        message += " xrms={:0.2f}".format(xrms)
        message += " yrms={:0.2f}".format(yrms)
        message += " xavg={:0.2f}".format(xavg)
        message += " yavg={:0.2f}".format(yavg)
        print(message)


        history["xrms"].append(xrms)
        history["yrms"].append(yrms)
        history["xavg"].append(xavg)
        history["yavg"].append(yavg)

    histories["bunch"] = copy.deepcopy(history)

    # Analysis
    # ------------------------------------------------------------------------------

    for history in histories.values():
        for key in history:
            history[key] = np.array(history[key])

    # Print errors
    for key in histories["envelope"]:
        deltas = histories["bunch"][key] - histories["envelope"][key]
        print("key:", key)
        print("max_abs_delta:", np.max(np.abs(deltas)))
        print("avg_abs_delta:", np.mean(np.abs(deltas)))

    # Plot rms bunch sizes
    for key in ["xrms", "yrms"]:
        fig, ax = plt.subplots(figsize=(5, 3))
        for i, model in enumerate(["envelope", "bunch"]):
            color = ["black", "red"][i]
            lw = [None, 0][i]
            ax.plot(histories[model][key], marker=".", lw=lw, color=color, label=model)
        ax.set_ylim(0.0, ax.get_ylim()[1] * 2.0)
        ax.set_xlabel("Turn")
        ax.set_ylabel("RMS [mm]")
        ax.legend(loc="upper right")
        plt.savefig(os.path.join(output_dir, f"fig_{key}"))
        plt.close()

    # Plot centroids
    for key in ["xavg", "yavg"]:
        fig, ax = plt.subplots(figsize=(5, 3))
        for i, model in enumerate(["envelope", "bunch"]):
            color = ["black", "red"][i]
            lw = [None, 0][i]
            ax.plot(histories[model][key], marker=".", lw=lw, color=color, label=model)
        ax.set_ylim(-5.0, 5.0)
        ax.set_xlabel("Turn")
        ax.set_ylabel("AVG [mm]")
        ax.legend(loc="upper right")
        plt.savefig(os.path.join(output_dir, f"fig_{key}"))
        plt.close()

    # Collect bunch/envelope data on final turn.
    particles = collect_bunch(bunch)["coords"]
    particles[:, :4] *= 1000.0

    env_cov_matrix = envelope.cov()
    env_cov_matrix[:4, :4] *= 1000.0**2

    env_centroid = envelope.centroid()
    env_centroid[:4] *= 1000.0

    xmax = 4.0 * np.std(particles, axis=0)
    limits = list(zip(-xmax, xmax))
    labels = ["x [mm]", "xp [mrad]", "y [mm]", "yp [mrad]", "z [m]", "dE [GeV]"]

    # Plot x-x'
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.hist2d(particles[:, 0], particles[:, 1], bins=100, range=[limits[0], limits[1]])
    plot_rms_ellipse(
        env_cov_matrix[0:2, 0:2],
        center=(env_centroid[0], env_centroid[1]),
        level=2.0,
        color="red",
        ax=ax,
    )
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    plt.savefig(os.path.join(output_dir, "fig_dist_x_xp"))
    plt.close()

    # Plot corner
    fig, axs = plot_corner(
        particles,
        limits=limits,
        bins=100,
        labels=labels,
    )
    for i in range(6):
        for j in range(i):
            env_cov_matrix_proj = project_cov_matrix(env_cov_matrix, axis=(j, i))
            plot_rms_ellipse(
                env_cov_matrix_proj,
                center=(env_centroid[j], env_centroid[i]),
                level=2.0,
                color="red",
                ax=axs[i, j],
            )
    plt.savefig(os.path.join(output_dir, "fig_dist_corner"))
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bunch-length", type=float, default=120.0)
    parser.add_argument("--kin-energy", type=float, default=1.300)
    parser.add_argument("--intensity", type=float, default=2e14)

    parser.add_argument("--dist", type=str, default="kv", choices=["kv", "waterbag", "gauss"])
    parser.add_argument("--mismatch-x", type=float, default=0.0)
    parser.add_argument("--mismatch-y", type=float, default=0.0)
    parser.add_argument("--offset-x", type=float, default=0.0)
    parser.add_argument("--offset-y", type=float, default=0.0)
    parser.add_argument("--tilt", type=float, default=0)

    parser.add_argument("--nparts", type=int, default=100_000)
    parser.add_argument("--turns", type=int, default=100)
    parser.add_argument("--sol", type=int, default=0)
    parser.add_argument("--sc", type=int, default=0)
    parser.add_argument("--sc-grid", type=int, default=64)

    parser.add_argument("--handle-unknown", type=str, default=None, choices=["drift", "fit"])
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(args)
