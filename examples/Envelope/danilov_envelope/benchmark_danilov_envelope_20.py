"""Benchmark {2, 2} Danilov envelope solver vs. PIC."""

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
from orbit.danilov_envelope import DanilovEnvelope20
from orbit.danilov_envelope import DanilovEnvelopeTracker20
from orbit.danilov_envelope import DanilovEnvelopeSolverNode20
from orbit.danilov_envelope import add_danilov_envelope_solver_nodes_20
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.space_charge.sc2p5d import SC2p5D_AccNode
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_fodo_lattice
from utils import BunchMonitor


# Parse arguments
# --------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--phase-adv-x", type=float, default=85.0)
parser.add_argument("--phase-adv-y", type=float, default=85.0)
parser.add_argument("--intensity", type=float, default=50.0)
parser.add_argument("--eps_x", type=float, default=10.00e-06)
parser.add_argument("--eps_y", type=float, default=10.00e-06)
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

envelope = DanilovEnvelope20(
    eps_x=args.eps_x,
    eps_y=args.eps_y,
    mass=mass_proton,
    kin_energy=1.000,
    length=100.0,
    intensity=(args.intensity * 1.00e14),
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

tracker = DanilovEnvelopeTracker20(lattice, path_length_max=args.max_part_length)
tracker.match_zero_sc(envelope)
envelope_init = envelope.copy()


# Track envelope
# --------------------------------------------------------------------------------------

histories = {}

envelope = envelope_init.copy()
history = tracker.track(envelope, history=True, periods=args.periods)
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

sc_calc = SpaceChargeCalc2p5D(64, 64, 1)
sc_path_length_min = 1.00e-06
sc_nodes = setSC2p5DAccNodes(lattice, sc_path_length_min, sc_calc)

bunch = envelope_init.to_bunch(env=False, size=100_000)
for periods in range(args.periods):
    lattice.trackBunch(bunch, actionContainer=action_container)

history = monitor.package_history()
histories["bunch"] = copy.deepcopy(history)


# Plot comparison
# --------------------------------------------------------------------------------------

figwidth = 4.0 * args.periods
figwidth = min(figwidth, 10.0)

fig, ax = plt.subplots(figsize=(figwidth, 2.5), constrained_layout=True)
for i, key in enumerate(histories):
    history = histories[key]

    plot_kws = {}
    if key == "envelope":
        plot_kws["ls"] = "-"
        plot_kws["lw"] = 2.5
        plot_kws["marker"] = None
    if key == "bunch":
        plot_kws["ls"] = "-"
        plot_kws["lw"] = 0.0
        plot_kws["marker"] = "."
        plot_kws["ms"] = 3.0
        plot_kws["color"] = "black"

    ax.plot(history["s"], history["xrms"] * 1000.0, **plot_kws)
    ax.plot(history["s"], history["yrms"] * 1000.0, **plot_kws)

ax.set_ylim(0.0, ax.get_ylim()[1])
ax.set_xlabel("Distance [m]")
ax.set_ylabel("Size [mm]")

filename = "fig_benchmark_rms.png"
filename = os.path.join(output_dir, filename)
plt.savefig(filename, dpi=300)
plt.show()
