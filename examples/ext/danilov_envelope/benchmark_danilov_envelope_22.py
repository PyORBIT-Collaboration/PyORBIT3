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
from orbit.danilov_envelope import DanilovEnvelope22
from orbit.danilov_envelope import DanilovEnvelopeTracker22
from orbit.danilov_envelope import DanilovEnvelopeSolverNode22
from orbit.danilov_envelope import add_danilov_envelope_solver_nodes_22
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
parser.add_argument("--intensity", type=float, default=20.0)
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

envelope = DanilovEnvelope22(
    intrinsic_emittance=20.00e-06,
    eps_x_frac=0.5,
    mass=0.938,
    kin_energy=1.0,
    length=100.0,
    intensity=args.intensity * 1.0e+14,
    mode=1,
    params=None
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

tracker = DanilovEnvelopeTracker22(lattice, path_length_max=args.max_part_length)
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

fig, axs = plt.subplots(figsize=(figwidth, 4.0), nrows=2, sharex=True, constrained_layout=True)
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

    axs[0].plot(history["s"], history["xrms"] * 1.00e+03, **plot_kws)
    axs[0].plot(history["s"], history["yrms"] * 1.00e+03, **plot_kws)
    axs[1].plot(history["s"], history["epsx"] * 1.00e+06, **plot_kws)
    axs[1].plot(history["s"], history["epsy"] * 1.00e+06, **plot_kws)

for ax in axs:
    ax.set_ylim(0.0, ax.get_ylim()[1])

axs[0].set_xlabel("Distance [m]")
axs[0].set_ylabel("Size [mm]")
axs[1].set_ylabel("Emittance [mm mrad]")

filename = "fig_benchmark.png"
filename = os.path.join(output_dir, filename)
plt.savefig(filename, dpi=300)
plt.show()
