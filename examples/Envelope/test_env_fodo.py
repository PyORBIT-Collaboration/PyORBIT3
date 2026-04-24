"""Test envelope tracker in FODO lattice."""
import argparse
import copy
import math

import numpy as np
import matplotlib.pyplot as plt

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.envelope import Envelope
from orbit.envelope import EnvelopeTracker
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

plt.style.use("style.mplstyle")


# Parse arguments
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--nslice", type=int, default=10)
parser.add_argument("--mismatch", type=float, default=0.0)
parser.add_argument("--turns", type=int, default=25)
parser.add_argument("--nparts", type=int, default=10_000)
args = parser.parse_args()


# Create lattice
# ------------------------------------------------------------------------------

nodes = [
    QuadTEAPOT(length=0.5, kq=+0.25),
    DriftTEAPOT(length=1.0),
    QuadTEAPOT(length=1.0, kq=-0.25),
    DriftTEAPOT(length=1.0),
    QuadTEAPOT(length=0.5, kq=+0.25),
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
sync_part.kinEnergy(1.0)

# Find periodic lattice parameters
matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
matrix_lattice_params = matrix_lattice.getRingParametersDict()
alpha_x = matrix_lattice_params["alpha x"]
alpha_y = matrix_lattice_params["alpha y"]
beta_x = matrix_lattice_params["beta x [m]"]
beta_y = matrix_lattice_params["beta y [m]"]
eps_x = 10.0e-06
eps_y = 10.0e-06

# Generate covariance matrix
cov_matrix = np.zeros((6, 6))
cov_matrix[0, 0] = eps_x * beta_x
cov_matrix[2, 2] = eps_y * beta_y
cov_matrix[0, 1] = cov_matrix[1, 0] = -eps_x * alpha_x
cov_matrix[2, 3] = cov_matrix[3, 2] = -eps_y * alpha_y
cov_matrix[1, 1] = eps_x * (1.0 + alpha_x**2) / beta_x
cov_matrix[3, 3] = eps_y * (1.0 + alpha_y**2) / beta_y
cov_matrix[4, 4] = 10.0 ** 2
cov_matrix[5, 5] = 0.0

# Mismatch x
cov_matrix[0, 0] *= (1.0 + args.mismatch) ** 2

# Create envelope
envelope = Envelope(
    sync_part=sync_part,
    cov_matrix=cov_matrix,
    centroid=np.zeros(6),
)


# Track envelope
# ------------------------------------------------------------------------------

print("TRACK ENVELOPE")

tracker = EnvelopeTracker(lattice)

history = {"xrms": [], "yrms": []}
for turn in range(args.turns):
    if turn > 0:
        tracker.track(envelope)

    cov_matrix = envelope.cov()
    xrms = 1000.0 * math.sqrt(cov_matrix[0, 0])
    yrms = 1000.0 * math.sqrt(cov_matrix[2, 2])
    print(f"turn={turn + 1} xrms={xrms:0.5f} yrms={yrms:0.5f}")

    history["xrms"].append(xrms)
    history["yrms"].append(yrms)

histories = {}
histories["envelope"] = copy.deepcopy(history)


# Track bunch
# ------------------------------------------------------------------------------

print("TRACK BUNCH")

rng = np.random.default_rng()
bunch_coords = rng.multivariate_normal(
    mean=np.zeros(6),
    cov=cov_matrix,
    size=args.nparts,
)
for i in range(bunch_coords.shape[0]):
    bunch.addParticle(*bunch_coords[i])

history = {"xrms": [], "yrms": []}
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
    print(f"turn={turn + 1} xrms={xrms:0.5f} yrms={yrms:0.5f}")

    history["xrms"].append(xrms)
    history["yrms"].append(yrms)

histories["bunch"] = copy.deepcopy(history)

# Analysis
# ------------------------------------------------------------------------------

for history in histories.values():
    for key in history:
        history[key] = np.array(history[key])

for key in histories["envelope"]:
    deltas = histories["bunch"][key] - histories["envelope"][key]
    print("key:", key)
    print("max_abs_delta:", np.max(np.abs(deltas)))
    print("avg_abs_delta:", np.mean(np.abs(deltas)))


# fig, ax = plt.subplots()