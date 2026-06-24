"""Test 3D envelope tracker in drift

The initial beam is a uniform-density ellipsoid in the x-y-z plane, with zero initial velocity.
The ellipsoid can have arbitrary size and orientation.
"""

import argparse
import copy
import math
import os
import pathlib

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.spacecharge import SpaceChargeCalc3D
from orbit.bunch_utils import collect_bunch
from orbit.envelope import Envelope
from orbit.envelope import EnvelopeTracker
from orbit.space_charge.sc3d import setSC3DAccNodes
from orbit.teapot import DriftTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.utils.consts import mass_proton

from plot import plot_rms_ellipse
from plot import plot_corner
from utils import gen_dist
from utils import project_cov_matrix

plt.style.use("style.mplstyle")


def rotation_matrix_3d(angle_x: float, angle_y: float, angle_z: float) -> np.ndarray:
    return scipy.spatial.transform.Rotation.from_euler("xyz", [angle_x, angle_y, angle_z]).as_matrix()


def build_cov_matrix_xyz(rms_sizes: np.ndarray, rotation_matrix: np.ndarray = None) -> np.ndarray:
    cov_matrix = np.diag(np.square(rms_sizes))
    if rotation_matrix is None:
        return cov_matrix
    return rotation_matrix @ cov_matrix @ rotation_matrix.T


def main(args: argparse.Namespace) -> None:

    # Setup
    # ------------------------------------------------------------------------------
    path = pathlib.Path(__file__)
    output_dir = os.path.join("outputs", path.stem)
    os.makedirs(output_dir, exist_ok=True)


    # Create lattice
    # ------------------------------------------------------------------------------
    node = DriftTEAPOT(length=args.length)
    node.setLength(args.length)
    node.setnParts(args.nslice)

    lattice = TEAPOT_Lattice()
    lattice.addNode(node)
    lattice.initialize()


    # Create envelope
    # ------------------------------------------------------------------------------
    bunch = Bunch()
    bunch.mass(mass_proton)
    sync_part = bunch.getSyncParticle()
    sync_part.kinEnergy(args.kin_energy)

    cov_matrix_init = np.zeros((6, 6))

    rotation_matrix = rotation_matrix_3d(
        math.radians(args.rot_x),
        math.radians(args.rot_y),
        math.radians(args.rot_z)
    )
    print(rotation_matrix)

    cov_matrix_xyz = build_cov_matrix_xyz([args.rms_x, args.rms_y, args.rms_z], rotation_matrix=rotation_matrix)

    lorentz_matrix = np.diag([1.0, 1.0, 1.0 / sync_part.gamma()])
    cov_matrix_xyz = lorentz_matrix @ cov_matrix_xyz @ lorentz_matrix.T

    spatial_idx = np.ix_([0, 2, 4], [0, 2, 4])
    cov_matrix_init[spatial_idx] = cov_matrix_xyz

    print(cov_matrix_xyz * 1e6)
    print()
    print(cov_matrix_init * 1e6)


    centroid_init = np.zeros(6)

    envelope = Envelope(
        bunch=bunch,
        cov_matrix=cov_matrix_init,
        centroid=centroid_init,
        intensity=args.intensity,
    )

    # Track envelope
    # ------------------------------------------------------------------------------
    print("TRACK ENVELOPE")

    tracker = EnvelopeTracker(lattice, space_charge=("3d" if args.sc else None))

    history = {"xrms": [], "yrms": [], "zrms": []}
    for turn in range(args.turns):
        if turn > 0:
            tracker.track(envelope)

        cov_matrix = envelope.cov()
        centroid = envelope.centroid()

        xrms = 1000.0 * math.sqrt(cov_matrix[0, 0])
        yrms = 1000.0 * math.sqrt(cov_matrix[2, 2])
        zrms = 1000.0 * math.sqrt(cov_matrix[4, 4]) * envelope.gamma()

        history["xrms"].append(xrms)
        history["yrms"].append(yrms)
        history["zrms"].append(zrms)

        print(f"turn={turn} xrms={xrms:0.3f} yrms={yrms:0.3f} zrms={zrms:0.3f}")

    histories = {}
    histories["envelope"] = copy.deepcopy(history)


    # Track bunch
    # ------------------------------------------------------------------------------
    print("TRACK BUNCH")

    bunch_coords = np.zeros((args.nparts, 6))
    bunch_coords[:, (0, 2, 4)] = gen_dist(
        size=args.nparts, cov_matrix=cov_matrix_xyz, name="waterbag"
    )

    for x, xp, y, yp, z, dE in bunch_coords:
        bunch.addParticle(x, xp, y, yp, z, dE)

    size_global = bunch.getSizeGlobal()
    bunch.macroSize(args.intensity / size_global)

    if args.sc:
        sc_calc = SpaceChargeCalc3D(args.sc_grid, args.sc_grid, args.sc_grid)
        sc_path_length_min = 0.01
        sc_nodes = setSC3DAccNodes(lattice, sc_path_length_min, sc_calc)

    history = {"xrms": [], "yrms": [], "zrms": []}
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
        zrms = 1000.0 * math.sqrt(cov_matrix[4, 4]) * bunch.getSyncParticle().gamma()

        history["xrms"].append(xrms)
        history["yrms"].append(yrms)
        history["zrms"].append(zrms)

        print(f"turn={turn} xrms={xrms:0.3f} yrms={yrms:0.3f} zrms={zrms:0.3f}")

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
    for key in ["xrms", "yrms", "zrms"]:
        fig, ax = plt.subplots(figsize=(4, 3))
        for i, model in enumerate(["envelope", "bunch"]):
            color = ["black", "red"][i]
            lw = [None, 0][i]
            ax.plot(histories[model][key], marker=".", lw=lw, color=color, label=model)
        ax.set_ylim(0.0, ax.get_ylim()[1] * 1.2)
        ax.set_xlabel("s [mm]")
        ax.set_ylabel("rms size [mm]")
        ax.legend(loc="upper left")
        plt.savefig(os.path.join(output_dir, f"fig_{key}"))
        plt.close()

    # Collect bunch/envelope data on final turn.
    particles = collect_bunch(bunch)["coords"]
    particles *= 1e3

    env_cov_matrix = envelope.cov()
    env_cov_matrix *= 1e6

    env_centroid = envelope.centroid()
    env_centroid *= 1e3

    xmax = 4.0 * np.std(particles, axis=0)
    limits = list(zip(-xmax, xmax))
    labels = ["x [mm]", "xp [mrad]", "y [mm]", "yp [mrad]", "z [mm]", "dE [MeV]"]

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kin-energy", type=float, default=0.0025)
    parser.add_argument("--intensity", type=float, default=3e10)

    parser.add_argument("--rms-x", type=float, default=0.010)
    parser.add_argument("--rms-y", type=float, default=0.010)
    parser.add_argument("--rms-z", type=float, default=0.010)

    parser.add_argument("--rot-x", type=float, default=0.0)
    parser.add_argument("--rot-y", type=float, default=0.0)
    parser.add_argument("--rot-z", type=float, default=0.0)

    parser.add_argument("--nslice", type=int, default=10)
    parser.add_argument("--length", type=float, default=0.1)
    parser.add_argument("--turns", type=int, default=20)
    parser.add_argument("--sc-grid", type=int, default=64)

    parser.add_argument("--nparts", type=int, default=100_000)
    parser.add_argument("--sc", type=int, default=0)
    args = parser.parse_args()

    main(args)