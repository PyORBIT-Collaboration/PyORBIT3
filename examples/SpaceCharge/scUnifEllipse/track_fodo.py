"""Test 2D envelope tracker in FODO lattice."""

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
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse2D
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import KVDist2D
from orbit.bunch_generators import WaterBagDist2D
from orbit.bunch_generators import GaussDist2D
from orbit.bunch_utils import collect_bunch
from orbit.lattice import AccNode
from orbit.lattice import AccLattice
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.space_charge.sc2p5d import setSCUnifEllipse2DAccNodes
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from plot import plot_rms_ellipse
from plot import plot_corner
from utils import build_rotation_matrix_xy
from utils import project_cov_matrix

plt.style.use("style.mplstyle")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bunch-length", type=float, default=5.0)
    parser.add_argument("--kin-energy", type=float, default=0.0025)
    parser.add_argument("--intensity", type=float, default=3e9)

    parser.add_argument(
        "--dist", type=str, default="kv", choices=["kv", "waterbag", "gauss"]
    )
    parser.add_argument("--mismatch-x", type=float, default=0.0)
    parser.add_argument("--mismatch-y", type=float, default=0.0)
    parser.add_argument("--tilt", type=float, default=0)

    parser.add_argument("--nslice", type=int, default=10)
    parser.add_argument("--kq", type=float, default=0.25)

    parser.add_argument("--nparts", type=int, default=100_000)
    parser.add_argument("--nturns", type=int, default=25)
    parser.add_argument("--sc", type=int, default=0)
    parser.add_argument("--sc-grid-res", type=int, default=128)
    parser.add_argument("--sc-ellipse-n", type=int, default=1)

    parser.add_argument("--plot-bins", type=int, default=64)
    parser.add_argument("--plot-mask", action="store_true")
    return parser.parse_args()


def make_lattice(args: argparse.Namespace) -> AccLattice:
    nodes = [
        QuadTEAPOT(length=0.5, kq=+args.kq),
        DriftTEAPOT(length=1.0),
        QuadTEAPOT(length=1.0, kq=-args.kq),
        DriftTEAPOT(length=1.0),
        QuadTEAPOT(length=0.5, kq=+args.kq),
    ]

    lattice = TEAPOT_Lattice()
    for node in nodes:
        node.setnParts(args.nslice)
        node.setUsageFringeFieldIN(False)
        node.setUsageFringeFieldOUT(False)
        lattice.addNode(node)
    lattice.initialize()
    return lattice


def make_empy_bunch(args: argparse.Namespace) -> Bunch:
    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(args.kin_energy)
    return bunch


def make_bunch(args: argparse.Namespace) -> Bunch:
    bunch = make_empy_bunch(args)
    lattice = make_lattice(args)

    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
    matrix_lattice_params = matrix_lattice.getRingParametersDict()

    alpha_x = matrix_lattice_params["alpha x"]
    alpha_y = matrix_lattice_params["alpha y"]
    beta_x = matrix_lattice_params["beta x [m]"]
    beta_y = matrix_lattice_params["beta y [m]"]
    eps_x = 0.25e-06
    eps_y = eps_x

    twiss_x = TwissContainer(alpha_x, beta_x, eps_x)
    twiss_y = TwissContainer(alpha_y, beta_y, eps_y)

    if args.dist == "kv":
        dist = KVDist2D(twiss_x, twiss_y)
    elif args.dist == "waterbag":
        dist = WaterBagDist2D(twiss_x, twiss_y)
    elif args.dist == "gauss":
        dist = GaussDist2D(twiss_x, twiss_y)
    else:
        raise ValueError

    particles = np.zeros((args.nparts, 6))
    for i in range(args.nparts):
        particles[i, :4] = dist.getCoordinates()
        particles[i, 4] = args.bunch_length * np.random.uniform(-1.0, 1.0)

    matrix = build_rotation_matrix_xy(math.radians(args.tilt))
    particles[:, :4] = particles[:, :4] @ matrix.T

    particles[:, 0] *= 1.0 + args.mismatch_x
    particles[:, 2] *= 1.0 + args.mismatch_y

    for i in range(args.nparts):
        bunch.addParticle(*particles[i])

    bunch_size = bunch.getSizeGlobal()
    bunch.macroSize(args.intensity / bunch_size)
    return bunch


def get_bunch_cov(bunch: Bunch) -> np.ndarray:
    twiss_calc = BunchTwissAnalysis()
    twiss_calc.analyzeBunch(bunch)

    cov_matrix = np.zeros((6, 6))
    for i in range(6):
        for j in range(i + 1):
            cov_matrix[i, j] = twiss_calc.getCorrelation(j, i)
            cov_matrix[j, i] = cov_matrix[i, j]
    return cov_matrix


def track(lattice: TEAPOT_Lattice, bunch: Bunch, nturns: int) -> dict:
    bunch_out = Bunch()
    bunch.copyBunchTo(bunch_out)

    history = {"rms_x": [], "rms_y": [], "eps_x": [], "eps_y": []}
    for turn in range(nturns):
        if turn > 0:
            lattice.trackBunch(bunch_out)

        cov_matrix = 1e6 * get_bunch_cov(bunch_out)
        x_rms = math.sqrt(cov_matrix[0, 0])
        y_rms = math.sqrt(cov_matrix[2, 2])
        eps_x = np.sqrt(np.linalg.det(cov_matrix[0:2, 0:2]))
        eps_y = np.sqrt(np.linalg.det(cov_matrix[2:4, 2:4]))

        history["rms_x"].append(x_rms)
        history["rms_y"].append(y_rms)
        history["eps_x"].append(eps_x)
        history["eps_y"].append(eps_y)

        print(
            f"turn={turn} xrms={x_rms:0.3f} yrms={y_rms:0.3f} epsx={eps_x:0.3f} epsy={eps_y:0.3f}"
        )

    particles_out = collect_bunch(bunch_out)["coords"]
    particles_out = particles_out[:, (0, 1, 2, 3)]
    particles_out *= 1000.0
    return {
        "particles": particles_out.copy(),
        "history": history,
    }


def main(args: argparse.Namespace) -> None:
    path = pathlib.Path(__file__)
    output_dir = os.path.join("outputs", path.stem, time.strftime("%Y%m%d_%H%M%S"))
    os.makedirs(output_dir, exist_ok=True)

    results = {}
    for key in ["SC2p5D", "SCUnifEllipse"]:
        results[key] = {}

    # Make initial bunch
    bunch = make_bunch(args)

    # Track bunch with SC2p5D space charge nodes
    lattice = make_lattice(args)
    if args.sc:
        sc_calc = SpaceChargeCalc2p5D(args.sc_grid_res, args.sc_grid_res, 1)
        sc_path_length_min = 1.00e-06
        sc_nodes = setSC2p5DAccNodes(lattice, sc_path_length_min, sc_calc)

    print("TRACK SC2p5D")
    results["SC2p5D"] = track(lattice, bunch, nturns=args.nturns)

    # Track bunch with SCUnifEllipse space charge nodes
    lattice = make_lattice(args)
    if args.sc:
        sc_calc = SpaceChargeCalcUnifEllipse2D(args.sc_ellipse_n)
        sc_nodes = setSCUnifEllipse2DAccNodes(lattice, sc_path_length_min, sc_calc)

    print("TRACK SCUnifEllipse2D")
    results["SCUnifEllipse"] = track(lattice, bunch, nturns=args.nturns)

    # Analysis
    # ------------------------------------------------------------------------------

    for model, result in results.items():
        history = result["history"]
        for key in history:
            history[key] = np.array(history[key])

    # Print errors
    for key in results["SCUnifEllipse"]["history"]:
        deltas = (
            results["SCUnifEllipse"]["history"][key] - results["SC2p5D"]["history"][key]
        )
        print("key:", key)
        print("max_abs_delta:", np.max(np.abs(deltas)))
        print("avg_abs_delta:", np.mean(np.abs(deltas)))

    # Plot rms size and emittance
    for key in ["eps_x", "eps_y", "rms_x", "rms_y"]:
        fig, ax = plt.subplots(figsize=(5, 3))
        for i, model in enumerate(["SC2p5D", "SCUnifEllipse"]):
            plot_kws = {}
            plot_kws["color"] = ["black", "red"][i]
            plot_kws["lw"] = [None, 0][i]
            ax.plot(results[model]["history"][key], marker=".", label=model, **plot_kws)
        ax.set_ylim(0.0, ax.get_ylim()[1] * 2.0)
        ax.set_xlabel("Turn")
        ax.set_ylabel(key)
        ax.legend(loc="upper right")
        plt.savefig(os.path.join(output_dir, f"fig_{key}"))
        plt.close()

    # Set plot limits
    particles = results["SC2p5D"]["particles"]
    xmax = 4.0 * np.std(particles, axis=0)
    limits = list(zip(-xmax, xmax))
    dims = ["x", "xp", "y", "yp"]
    labels = ["x [mm]", "xp [mrad]", "y [mm]", "yp [mrad]"]

    # Plot x-x', y-y', x-y
    for axis in [(0, 1), (2, 3), (0, 2)]:
        fig, axs = plt.subplots(figsize=(6, 3), ncols=2, sharex=True, sharey=True)
        for ax, model in zip(axs, results):
            particles = results[model]["particles"]
            values, edges = np.histogramdd(
                particles[:, axis], bins=args.plot_bins, range=[limits[k] for k in axis]
            )
            if args.plot_mask:
                values = np.ma.masked_equal(values, 0.0)
            ax.pcolormesh(edges[0], edges[1], values.T)
            ax.set_xlabel(labels[axis[0]])
            ax.set_ylabel(labels[axis[1]])
            ax.set_title(model)
        plt.savefig(
            os.path.join(output_dir, f"fig_dist_{dims[axis[0]]}_{dims[axis[1]]}")
        )
        plt.close()

    # Plot corner
    for model in results:
        particles = results[model]["particles"]
        fig, axs = plot_corner(
            particles,
            limits=limits,
            bins=args.plot_bins,
            labels=labels,
            mask=args.plot_mask,
        )
        plt.savefig(os.path.join(output_dir, f"fig_dist_corner_{model}"))
        plt.close()


if __name__ == "__main__":
    main(parse_args())
