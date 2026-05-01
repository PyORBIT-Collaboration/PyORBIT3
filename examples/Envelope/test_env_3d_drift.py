"""Test 3D envelope tracker in drift."""

import argparse
import copy
import math
import os
import pathlib

import numpy as np
import matplotlib.pyplot as plt

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


# Parse arguments
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--kin-energy", type=float, default=0.025)
parser.add_argument("--intensity", type=float, default=5e11)

parser.add_argument("--xrms", type=float, default=0.010)
parser.add_argument("--yrms", type=float, default=0.010)
parser.add_argument("--zrms", type=float, default=0.010)

parser.add_argument("--nslice", type=int, default=10)
parser.add_argument("--length", type=float, default=0.1)
parser.add_argument("--turns", type=int, default=20)
parser.add_argument("--sc-grid", type=int, default=64)

parser.add_argument("--nparts", type=int, default=100_000)
parser.add_argument("--sc", type=int, default=0)
args = parser.parse_args()


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
cov_matrix_init[0, 0] = args.xrms**2
cov_matrix_init[2, 2] = args.yrms**2
cov_matrix_init[4, 4] = (args.zrms / sync_part.gamma()) ** 2

centroid_init = np.zeros(6)

envelope = Envelope(
    sync_part=sync_part,
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
    args.nparts, cov_matrix=np.eye(3), name="waterbag"
)
bunch_coords[:, 0] *= args.xrms
bunch_coords[:, 2] *= args.yrms
bunch_coords[:, 4] *= args.zrms / sync_part.gamma()

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
    fig, ax = plt.subplots(figsize=(5, 3))
    for i, model in enumerate(["envelope", "bunch"]):
        color = ["black", "red"][i]
        lw = [None, 0][i]
        ax.plot(histories[model][key], marker=".", lw=lw, color=color, label=model)
    ax.set_ylim(0.0, ax.get_ylim()[1])
    ax.set_xlabel("Turn")
    ax.set_ylabel("RMS [mm]")
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
