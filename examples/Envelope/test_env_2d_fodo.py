"""Test 2D envelope tracker in FODO lattice."""

import argparse
import copy
import math
import os
import pathlib

import numpy as np
import matplotlib.pyplot as plt

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.bunch_utils import collect_bunch
from orbit.envelope import Envelope
from orbit.envelope import EnvelopeTracker
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from plot import plot_rms_ellipse
from plot import plot_corner
from utils import gen_dist
from utils import build_rotation_matrix_xy
from utils import project_cov_matrix

plt.style.use("style.mplstyle")


# Parse arguments
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--zrms", type=float, default=5.0)
parser.add_argument("--kin-energy", type=float, default=0.0025)
parser.add_argument("--intensity", type=float, default=5e9)

parser.add_argument(
    "--dist", type=str, default="kv", choices=["kv", "waterbag", "gauss"]
)
parser.add_argument("--mismatch-x", type=float, default=0.0)
parser.add_argument("--mismatch-y", type=float, default=0.0)
parser.add_argument("--offset-x", type=float, default=0.0)
parser.add_argument("--offset-y", type=float, default=0.0)
parser.add_argument("--tilt", type=float, default=0)

parser.add_argument("--nslice", type=int, default=10)
parser.add_argument("--kq", type=float, default=0.25)

parser.add_argument("--nparts", type=int, default=100_000)
parser.add_argument("--turns", type=int, default=25)
parser.add_argument("--sc", type=int, default=0)
args = parser.parse_args()


# Setup
# ------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


# Create lattice
# ------------------------------------------------------------------------------

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
eps_x = 0.25e-06
eps_y = eps_x

# Generate covariance matrix
cov_matrix = np.zeros((6, 6))
cov_matrix[0, 0] = eps_x * beta_x
cov_matrix[2, 2] = eps_y * beta_y
cov_matrix[0, 1] = cov_matrix[1, 0] = -eps_x * alpha_x
cov_matrix[2, 3] = cov_matrix[3, 2] = -eps_y * alpha_y
cov_matrix[1, 1] = eps_x * (1.0 + alpha_x**2) / beta_x
cov_matrix[3, 3] = eps_y * (1.0 + alpha_y**2) / beta_y
cov_matrix[4, 4] = args.zrms**2
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
    sync_part=sync_part,
    cov_matrix=cov_matrix_init,
    centroid=centroid_init,
    intensity=args.intensity,
)


# Track envelope
# ------------------------------------------------------------------------------

print("TRACK ENVELOPE")

tracker = EnvelopeTracker(lattice, space_charge=("2d" if args.sc else None))

history = {"xrms": [], "yrms": [], "xavg": [], "yavg": []}
for turn in range(args.turns):
    if turn > 0:
        tracker.track(envelope)

    cov_matrix = envelope.cov()
    centroid = envelope.centroid()

    xrms = 1000.0 * math.sqrt(cov_matrix[0, 0])
    yrms = 1000.0 * math.sqrt(cov_matrix[2, 2])
    xavg = 1000.0 * centroid[0]
    yavg = 1000.0 * centroid[2]

    print(
        f"turn={turn} xrms={xrms:0.3f} yrms={yrms:0.3f} xavg={xavg:0.3f} yavg={yavg:0.3f}"
    )

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
    n=args.nparts, cov_matrix=cov_matrix_init[0:4, 0:4], name=args.dist
)
bunch_coords[:, 4] = 2.0 * rng.uniform(-args.zrms, args.zrms, size=args.nparts)
bunch_coords += centroid_init[None, :6]

for i in range(bunch_coords.shape[0]):
    bunch.addParticle(*bunch_coords[i])

if args.sc:
    sc_calc = SpaceChargeCalc2p5D(128, 128, 1)
    sc_path_length_min = 1.00e-06
    sc_nodes = setSC2p5DAccNodes(lattice, sc_path_length_min, sc_calc)

    bunch_size = bunch.getSizeGlobal()
    bunch.macroSize(args.intensity / bunch_size)

history = {"xrms": [], "yrms": [], "xavg": [], "yavg": []}
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

    print(
        f"turn={turn} xrms={xrms:0.3f} yrms={yrms:0.3f} xavg={xavg:0.3f} yavg={yavg:0.3f}"
    )

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
