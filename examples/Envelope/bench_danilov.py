"""Benchmark Danilov envelope tracker vs. PIC."""

import argparse
import copy
import os
import pathlib
from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.envelope import DanilovEnvelope
from orbit.envelope import DanilovEnvelopeMonitor
from orbit.envelope import DanilovEnvelopeTracker
from orbit.envelope import add_danilov_envelope_tracker_nodes
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.space_charge.sc2p5d import SC2p5D_AccNode
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_fodo_lattice
from utils import BunchMonitor

plt.style.use("style.mplstyle")


# Parse arguments
# --------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--phase-adv-x", type=float, default=85.0)
parser.add_argument("--phase-adv-y", type=float, default=85.0)
parser.add_argument("--intensity", type=float, default=20.0e14)
parser.add_argument("--max-part-length", type=float, default=0.1)
parser.add_argument("--mismatch", type=float, default=0.0)
parser.add_argument("--periods", type=int, default=5)
args = parser.parse_args()


# Setup
# --------------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


# Set up simulation
# --------------------------------------------------------------------------------------

envelope = DanilovEnvelope(
    eps_1=20.00e-06,
    eps_2=0.0,
    mass=0.938,
    kin_energy=1.0,
    length=100.0,
    line_density=(args.intensity / 100.0),
    params=None,
)

lattice = make_fodo_lattice(
    phase_adv_x=np.radians(args.phase_adv_x),
    phase_adv_y=np.radians(args.phase_adv_y),
    length=5.0,
    mass=envelope.mass,
    kin_energy=envelope.kin_energy,
    max_part_length=args.max_part_length,
    verbose=1,
)

tracker = DanilovEnvelopeTracker(lattice, path_length_max=args.max_part_length)
tracker.match_zero_sc(envelope, method="2d")
envelope_init = envelope.copy()


# Track envelope
# --------------------------------------------------------------------------------------

histories = {}

history = tracker.track(envelope_init.copy(), periods=args.periods, history=True)
histories["envelope"] = copy.deepcopy(history)


# Track bunch
# --------------------------------------------------------------------------------------

lattice = make_fodo_lattice(
    phase_adv_x=np.radians(args.phase_adv_x),
    phase_adv_y=np.radians(args.phase_adv_y),
    length=5.0,
    mass=envelope.mass,
    kin_energy=envelope.kin_energy,
    max_part_length=args.max_part_length,
    verbose=1,
)

monitor = BunchMonitor()
action_container = AccActionsContainer()
action_container.addAction(monitor, AccActionsContainer.ENTRANCE)
action_container.addAction(monitor, AccActionsContainer.EXIT)

sc_calc = SpaceChargeCalc2p5D(128, 128, 1)
sc_path_length_min = 1.00e-06
sc_nodes = setSC2p5DAccNodes(lattice, sc_path_length_min, sc_calc)

bunch = envelope_init.to_bunch(env=False, size=128_000)

for periods in range(args.periods):
    lattice.trackBunch(bunch, actionContainer=action_container)

history = monitor.package_history()
histories["bunch"] = copy.deepcopy(history)


# Plot comparison
# --------------------------------------------------------------------------------------

figwidth = 3.0 * args.periods
figwidth = min(figwidth, 7.0)

fig, axs = plt.subplots(nrows=2, figsize=(figwidth, 4.0), sharex=True, sharey=True)
for i, ax in enumerate(axs):
    param = ["xrms", "yrms"][i]
    for j, key in enumerate(histories):
        history = histories[key]
        if key == "envelope":
            ax.plot(
                history["s"],
                np.multiply(history[param], 1000.0),
                color="black",
                lw=1.5,
            )
        else:
            stride = 10
            ax.plot(
                history["s"][::stride],
                np.multiply(history[param][::stride], 1000.0),
                marker=".",
                lw=0,
                color="red",
            )

for ax in axs:
    ax.set_ylim(0.0, ax.get_ylim()[1])
axs[1].set_xlabel("Distance [m]")
axs[0].set_ylabel("RMS x [mm]")
axs[1].set_ylabel("RMS y [mm]")

filename = "fig_benchmark_rms.png"
filename = os.path.join(output_dir, filename)
plt.savefig(filename, dpi=300)
plt.show()
