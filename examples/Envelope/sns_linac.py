"""SNS linac benchmark."""

import argparse
import math
import os
import random
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.linac import BaseRfGap
from orbit.core.linac import MatrixRfGap
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse
from orbit.core.spacecharge import SpaceChargeCalc3D
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D
from orbit.bunch_generators import GaussDist3D
from orbit.bunch_generators import KVDist3D
from orbit.bunch_utils import collect_bunch
from orbit.envelope import Envelope
from orbit.envelope import EnvelopeAnalysis
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory
from orbit.py_linac.lattice import LinacAccLattice
from orbit.space_charge.sc3d import setSC3DAccNodes
from orbit.space_charge.sc3d import setUniformEllipsesSCAccNodes
from orbit.utils.consts import mass_proton
from orbit.utils.consts import mass_electron
from orbit.utils.consts import charge_electron

from plot import plot_corner
from plot import plot_rms_ellipse
from utils import project_cov_matrix

plt.style.use("style.mplstyle")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dist", type=str, default="kv")
    parser.add_argument("--nparts", type=int, default=20_000)
    parser.add_argument("--sc", type=int, default=0)
    parser.add_argument(
        "--sc-model", type=str, default="ellipsoid", choices=("ellipsoid", "fft")
    )
    parser.add_argument("--sc-ellipsoid-n", type=int, default=1)
    parser.add_argument("--sc-path-length-min", type=float, default=0.01)
    parser.add_argument("--current", type=float, default=0.038)
    parser.add_argument("--show", type=int, default=0)
    parser.add_argument("--seq-stop", type=str, default="CCL2")
    parser.add_argument(
        "--fit", type=int, default=1, help="Use best-fit transfer matrices"
    )
    return parser.parse_args()


def make_lattice(args: argparse.Namespace) -> LinacAccLattice:
    seq_names = [
        "MEBT",
        "DTL1",
        "DTL2",
        "DTL3",
        "DTL4",
        "DTL5",
        "DTL6",
        "CCL1",
        "CCL2",
        "CCL3",
        "CCL4",
        "SCLMed",
        "SCLHigh",
        "HEBT1",
        "HEBT2",
    ]
    if args.seq_stop:
        index = seq_names.index(args.seq_stop) + 1
        seq_names = seq_names[:index]

    sns_linac_factory = SNS_LinacLatticeFactory()
    sns_linac_factory.setMaxDriftLength(args.sc_path_length_min)
    lattice = sns_linac_factory.getLinacAccLattice(seq_names, "inputs/sns_linac.xml")

    for node in lattice.getNodes():
        try:
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
        except:
            pass

    rf_gaps = lattice.getRF_Gaps()
    for rf_gap in rf_gaps:
        rf_gap.setCppGapModel(MatrixRfGap())

    return lattice


def make_bunch(args: argparse.Namespace) -> Bunch:
    kin_energy = 0.0025  # [GeV]
    mass = mass_proton + 2.0 * mass_electron
    frequency = 402.5e06
    charge = -1.0
    intensity = args.current / frequency / (math.fabs(charge) * charge_electron)

    bunch = Bunch()
    bunch.mass(mass)
    bunch.macroSize(intensity / args.nparts)
    bunch.charge(charge)

    sync_part = bunch.getSyncParticle()
    sync_part.kinEnergy(kin_energy)
    sync_part.time(0.0)

    alpha_x, beta_x, eps_x = (-1.962, 0.183, 2.874e-06)
    alpha_y, beta_y, eps_y = (+1.768, 0.162, 2.874e-06)
    alpha_z, beta_z, eps_z = (-0.0196, 116.414, 1.651e-08)

    twiss_x = TwissContainer(alpha_x, beta_x, eps_x)
    twiss_y = TwissContainer(alpha_y, beta_y, eps_y)
    twiss_z = TwissContainer(alpha_z, beta_z, eps_z)

    if args.dist == "waterbag":
        dist = WaterBagDist3D(twiss_x, twiss_y, twiss_z)
    elif args.dist == "kv":
        dist = KVDist3D(twiss_x, twiss_y, twiss_z)
    elif args.dist == "gauss":
        dist = GaussDist3D(twiss_x, twiss_y, twiss_z)
    else:
        raise ValueError("Unknown distribution '{}'".format(args.dist))

    for _ in range(args.nparts):
        bunch.addParticle(*dist.getCoordinates())
    return bunch


class EnvelopeMonitor:
    def __init__(self) -> None:
        self.history = {}

    def __call__(self, params_dict: dict) -> None:
        envelope = params_dict["envelope"]
        s = params_dict["path_length"]
        if self.history and self.history["s"][-1] == s:
            return
        parameters = {"s": s, **EnvelopeAnalysis(envelope).to_dict()}
        for key, value in parameters.items():
            self.history.setdefault(key, []).append(value)

        message = ""
        message += " s={:0.3f}".format(s)
        message += " rms_x={:0.3f}".format(1000.0 * parameters["rms_x"])
        message += " rms_y={:0.3f}".format(1000.0 * parameters["rms_y"])
        message += " rms_z={:0.3f}".format(1000.0 * parameters["rms_z"])
        message += " node={}".format(params_dict["node"].getName())
        print(message)


class BunchMonitor:
    def __init__(self) -> None:
        self.position_start = 0.0

        keys = [
            "s",
            "kin_energy",
            "gamma",
            "beta",
            "rms_x",
            "rms_y",
            "rms_z",
            "eps_x",
            "eps_x_n",
            "eps_y",
            "eps_y_n",
            "eps_z",
            "eps_z_n",
            "eps_1",
            "eps_2",
            "eps_3",
        ]
        self.history = {key: [] for key in keys}

    def __call__(self, params_dict: dict) -> None:
        bunch = params_dict["bunch"]
        node = params_dict["node"]
        position = params_dict["path_length"]

        if params_dict["old_pos"] == position:
            return
        if params_dict["old_pos"] + params_dict["pos_step"] > position:
            return
        params_dict["old_pos"] = position
        params_dict["count"] += 1

        analysis = EnvelopeAnalysis(Envelope(bunch=bunch))
        results = {"s": position, **analysis.to_dict()}
        for key in self.history:
            self.history[key].append(results[key])

        message = ""
        message += " s={:0.3f}".format(position + self.position_start)
        message += " rms_x={:0.3f}".format(1000.0 * results["rms_x"])
        message += " rms_y={:0.3f}".format(1000.0 * results["rms_y"])
        message += " rms_z={:0.3f}".format(1000.0 * results["rms_z"])
        message += " node={}".format(node.getName())
        print(message)


def make_output_dir(args: argparse.Namespace) -> str:
    output_dir = time.strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join("outputs", Path(__file__).stem, output_dir)
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def main(args: argparse.Namespace) -> None:
    output_dir = make_output_dir(args)
    random.seed(23)

    # Track envelope
    bunch = make_bunch(args)
    envelope = Envelope(bunch=bunch)

    lattice = make_lattice(args)
    lattice.trackDesignBunch(bunch)

    monitor = EnvelopeMonitor()
    action_container = AccActionsContainer()
    action_container.addAction(monitor, AccActionsContainer.ENTRANCE)
    action_container.addAction(monitor, AccActionsContainer.EXIT)
    lattice.trackEnvelope(
        envelope,
        action_container=action_container,
        sc=("3d" if args.sc else None),
        fit=args.fit,
    )
    histories = {"envelope": monitor.history}

    # Track bunch
    bunch = make_bunch(args)
    lattice = make_lattice(args)
    lattice.trackDesignBunch(bunch)
    if args.sc:
        if args.sc_model == "ellipsoid":
            n_ellipsoids = args.sc_ellipsoid_n
            sc_calc = SpaceChargeCalcUnifEllipse(n_ellipsoids)
            sc_nodes = setUniformEllipsesSCAccNodes(
                lattice, args.sc_path_length_min, sc_calc
            )
        if args.sc_model == "fft":
            sc_calc = SpaceChargeCalc3D(64, 64, 64)
            sc_nodes = setSC3DAccNodes(lattice, args.sc_path_length_min, sc_calc)

    monitor = BunchMonitor()
    action_container = AccActionsContainer()
    action_container.addAction(monitor, AccActionsContainer.ENTRANCE)
    action_container.addAction(monitor, AccActionsContainer.EXIT)
    params_dict = {"old_pos": -1.0, "count": 0, "pos_step": args.sc_path_length_min}
    lattice.trackBunch(bunch, paramsDict=params_dict, actionContainer=action_container)
    histories["bunch"] = monitor.history

    # Analysis
    # --------------------------------------------------------------------------------

    for mode in histories:
        for key in histories[mode]:
            histories[mode][key] = np.array(histories[mode][key])

    # History: rms size
    plot_kws = {}
    plot_kws["bunch"] = dict(color="black", lw=0, marker=".", ms=2)
    plot_kws["envelope"] = dict(color="red", lw=0, marker=".", ms=1)

    fig, axs = plt.subplots(
        nrows=3, figsize=(10, 5), sharex=True, constrained_layout=True
    )
    for mode in ["bunch", "envelope"]:
        history = histories[mode]
        for ax, key in zip(axs, ["rms_x", "rms_y", "rms_z"]):
            ax.plot(history["s"], history[key], **plot_kws[mode], label=mode)
    for ax in axs:
        ax.legend(loc="lower right")
    axs[0].set_ylabel("x rms [mm]")
    axs[1].set_ylabel("y rms [mm]")
    axs[2].set_ylabel("z rms [mm]")
    axs[2].set_xlabel("s [m]")
    plt.savefig(os.path.join(output_dir, "fig_history_rms_size.png"))
    if args.show:
        plt.show()
    plt.close()

    # History: rms emittance
    fig, axs = plt.subplots(
        nrows=3, figsize=(10, 5), sharex=True, constrained_layout=True
    )
    for mode in ["bunch", "envelope"]:
        history = histories[mode]
        for ax, key in zip(axs, ["eps_x_n", "eps_y_n", "eps_z_n"]):
            ax.plot(history["s"], history[key], **plot_kws[mode], label=mode)
    for ax in axs:
        ax.legend(loc="lower right")
        ax.set_ylim(0.0, ax.get_ylim()[1] * 2.0)
    axs[0].set_ylabel(r"$\gamma \beta \varepsilon_x$ [m rad]")
    axs[1].set_ylabel(r"$\gamma \beta  \varepsilon_y$ [m rad]")
    axs[2].set_ylabel(r"$\varepsilon_z / \beta$ [m GeV]")
    axs[2].set_xlabel("s [m]")
    plt.savefig(os.path.join(output_dir, "fig_history_rms_emittance.png"))
    if args.show:
        plt.show()
    plt.close()

    # History: energy
    fig, ax = plt.subplots(figsize=(5, 3))
    for mode in ["bunch", "envelope"]:
        history = histories[mode]
        ax.plot(history["s"], history["kin_energy"], **plot_kws[mode], label=mode)
    ax.legend(loc="lower right")
    ax.set_ylabel("energy [GeV]")
    ax.set_xlabel("s [m]")
    plt.savefig(os.path.join(output_dir, "fig_history_energy.png"))
    plt.close()

    # Collect bunch/envelope data
    particles = collect_bunch(bunch)["coords"]
    particles *= 1e3

    env_cov_matrix = EnvelopeAnalysis(envelope).covariance()
    env_cov_matrix *= 1e6

    env_centroid = envelope.centroid
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
    main(parse_args())
