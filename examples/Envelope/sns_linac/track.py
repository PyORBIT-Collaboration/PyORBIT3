import argparse
import math
import os
import random
import time

import numpy as np
import matplotlib.pyplot as plt

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.bunch import SyncParticle
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
from orbit.envelope import EnvelopeTracker
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

# local
from diagnostics import BunchMonitor

plt.style.use("style.mplstyle")


# Parse arguments
# --------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument(
    "--seq",
    type=str,
    default=None,
    choices=[
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
    ],
)
parser.add_argument("--sc", type=int, default=0)
parser.add_argument("--sc-model", type=str, default="ellipsoid")
parser.add_argument("--nparts", type=int, default=10_000)
parser.add_argument("--current", type=float, default=0.038)
parser.add_argument("--sc-path-length-min", type=float, default=0.01)
parser.add_argument("--show", type=int, default=0)
args = parser.parse_args()


# Setup
# --------------------------------------------------------------------------------

output_dir = "outputs"
os.makedirs(output_dir, exist_ok=True)

random.seed(100)


# Bunch
# --------------------------------------------------------------------------------

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

dist = WaterBagDist3D(twiss_x, twiss_y, twiss_z)
for _ in range(args.nparts):
    bunch.addParticle(*dist.getCoordinates())


# Lattice
# --------------------------------------------------------------------------------

sequence_names = [
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
if args.seq:
    stop_index = sequence_names.index(args.seq) + 1
    sequence_names = sequence_names[:stop_index]

sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)
lattice = sns_linac_factory.getLinacAccLattice(sequence_names, "sns_linac.xml")

for node in lattice.getNodes():
    try:
        node.setUsageFringeFieldIN(False)
        node.setUsageFringeFieldOUT(False)
    except:
        pass

rf_gaps = lattice.getRF_Gaps()
for rf_gap in rf_gaps:
    rf_gap.setCppGapModel(MatrixRfGap())

for index, node in enumerate(lattice.getNodes()):
    print(index, type(node), node.getName())

lattice.trackDesignBunch(bunch)


# Track envelope
# --------------------------------------------------------------------------------

twiss_calc = BunchTwissAnalysis()
twiss_calc.analyzeBunch(bunch)

cov_matrix = np.zeros((6, 6))
for i in range(6):
    for j in range(6):
        cov_matrix[i, j] = cov_matrix[j, i] = twiss_calc.getCorrelation(i, j)

envelope = Envelope(bunch=bunch, cov_matrix=cov_matrix, intensity=intensity)

envelope_tracker = EnvelopeTracker(lattice, space_charge=("3d" if args.sc else None))

histories = {}
histories["envelope"] = envelope_tracker.track_history(envelope)


# Track bunch
# --------------------------------------------------------------------------------

if args.sc:
    sc_path_length_min = args.sc_path_length_min
    if args.sc_model == "ellipsoid":
        n_ellipsoids = 1
        sc_calc = SpaceChargeCalcUnifEllipse(n_ellipsoids)
        sc_nodes = setUniformEllipsesSCAccNodes(lattice, sc_path_length_min, sc_calc)
    if args.sc_model == "3d":
        sc_calc = SpaceChargeCalc3D(64, 64, 64)
        sc_nodes = setSC3DAccNodes(lattice, sc_path_length_min, sc_calc)


monitor = BunchMonitor()

action_container = AccActionsContainer()
action_container.addAction(monitor, AccActionsContainer.ENTRANCE)
action_container.addAction(monitor, AccActionsContainer.EXIT)

params_dict = {"old_pos": -1.0, "count": 0, "pos_step": 0.1}

lattice.trackBunch(bunch, paramsDict=params_dict, actionContainer=action_container)

histories["bunch"] = monitor.history


# Analysis
# --------------------------------------------------------------------------------


# History: rms
for mode in histories:
    for key in histories[mode]:
        histories[mode][key] = np.array(histories[mode][key])

plot_kws = {}
plot_kws["bunch"] = dict(
    color="black",
    ls="-",
)
plot_kws["envelope"] = dict(
    color="red",
    # ls="--",
    lw=0,
    marker=".",
    ms=1,
)

fig, axs = plt.subplots(nrows=3, figsize=(5, 7), sharex=True, constrained_layout=True)
for mode in ["bunch", "envelope"]:
    history = histories[mode]
    for ax, key in zip(axs, ["rms_x", "rms_y", "rms_z"]):
        ax.plot(history["position"], history[key], **plot_kws[mode], label=mode)
for ax in axs:
    ax.legend(loc="lower right")
axs[0].set_ylabel("x rms [mm]")
axs[1].set_ylabel("y rms [mm]")
axs[2].set_ylabel("z rms [mm]")
axs[2].set_xlabel("s [m]")
plt.savefig(os.path.join(output_dir, "fig_history_rms.png"))
if args.show:
    plt.show()
plt.close()

# History: energy
fig, ax = plt.subplots(figsize=(5, 3))
for mode in ["bunch", "envelope"]:
    history = histories[mode]
    ax.plot(history["position"], history["kin_energy"], **plot_kws[mode], label=mode)
ax.legend(loc="lower right")
ax.set_ylabel("energy [GeV]")
ax.set_xlabel("s [m]")
plt.savefig(os.path.join(output_dir, "fig_history_energy.png"))
plt.close()