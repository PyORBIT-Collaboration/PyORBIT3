import argparse
import math
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
from orbit.utils.consts import charge_electron
from orbit.utils.consts import speed_of_light

plt.style.use("style.mplstyle")


# Parse arguments
# --------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--seq", type=str, default="MEBT")
parser.add_argument("--sc", type=int, default=0)
parser.add_argument("--sc-model", type=str, default="ellipsoid")
parser.add_argument("--nparts", type=int, default=10_000)
parser.add_argument("--current", type=float, default=0.038)
args = parser.parse_args()


# Setup
# --------------------------------------------------------------------------------

random.seed(100)


# Bunch
# --------------------------------------------------------------------------------

kin_energy = 0.0025  # [GeV]
mass = mass_proton
frequency = 402.5e+06
charge = -1.0
intensity = args.current / frequency / (math.fabs(charge) * charge_electron)

bunch = Bunch()
bunch.mass(mass)
bunch.getSyncParticle().kinEnergy(kin_energy)
bunch.macroSize(intensity / args.nparts)
bunch.charge(charge)

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

envelope = Envelope(bunch=bunch, cov_matrix=cov_matrix)

envelope_tracker = EnvelopeTracker(lattice, space_charge=None)

histories = {}
histories["envelope"] = envelope_tracker.track_history(envelope)


# Track bunch
# --------------------------------------------------------------------------------

lattice.trackDesignBunch(bunch)

if args.sc:
    sc_path_length_min = 0.01
    if args.sc_model == "ellipsoid":
        n_ellipsoids = 1
        sc_calc = SpaceChargeCalcUnifEllipse(n_ellipsoids)
        sc_nodes = setUniformEllipsesSCAccNodes(lattice, sc_path_length_min, sc_calc)
    if args.sc_model == "3d":
        sc_calc = SpaceChargeCalc3D(64, 64, 64)
        sc_nodes = setSC3DAccNodes(lattice, sc_path_length_min, sc_calc)


class BunchMonitor:
    def __init__(self) -> None:
        self.twiss_calc = BunchTwissAnalysis()
        self.position_start = 0.0

        self.history = {}
        self.history["position"] = []
        self.history["rms_x"] = []
        self.history["rms_y"] = []
        self.history["rms_z"] = []

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

        self.twiss_calc.analyzeBunch(bunch)

        cov_matrix = np.zeros((6, 6))
        for i in range(6):
            for j in range(6):
                cov_matrix[i, j] = cov_matrix[j, i] = self.twiss_calc.getCorrelation(i, j)

        xrms = 1000.0 * np.sqrt(cov_matrix[0, 0])
        yrms = 1000.0 * np.sqrt(cov_matrix[2, 2])
        zrms = 1000.0 * np.sqrt(cov_matrix[4, 4])

        message = ""
        message += " s={:0.3f}".format(position + self.position_start)
        message += " xrms={:0.3f}".format(xrms)
        message += " yrms={:0.3f}".format(yrms)
        message += " zrms={:0.3f}".format(zrms)
        message += " node={}".format(node.getName())
        print(message)

        self.history["position"].append(position + self.position_start)
        self.history["rms_x"].append(xrms)
        self.history["rms_y"].append(yrms)
        self.history["rms_z"].append(zrms)


monitor = BunchMonitor()

action_container = AccActionsContainer()
action_container.addAction(monitor, AccActionsContainer.ENTRANCE)
action_container.addAction(monitor, AccActionsContainer.EXIT)

params_dict = {"old_pos": -1.0, "count": 0, "pos_step": 0.1}

lattice.trackBunch(bunch, paramsDict=params_dict, actionContainer=action_container)

histories["bunch"] = monitor.history


# Analysis
# --------------------------------------------------------------------------------


# History
for mode in histories:
    for key in histories[mode]:
        histories[mode][key] = np.array(histories[mode][key])

fig, axs = plt.subplots(nrows=3, figsize=(5, 7), sharex=True, constrained_layout=True)
for i, mode in enumerate(["bunch", "envelope"]):
    history = histories[mode]
    color = ["black", "red"][i]
    ls = ["-", "--"][i]
    for ax, key in zip(axs, ["rms_x", "rms_y", "rms_z"]):
        ax.plot(history["position"], history[key], color=color, ls=ls, label=mode)

for ax in axs:
    ax.legend(loc="lower right")
axs[0].set_ylabel("x rms [mm]")
axs[1].set_ylabel("y rms [mm]")
axs[2].set_ylabel("z rms [mm]")
axs[2].set_xlabel("s [m]")
plt.show()


# Final coordinates
bunch_coords = collect_bunch(bunch)["coords"]
bunch_cov_matrix = np.cov(bunch_coords.T)

print(np.round(1000.0 * np.sqrt(np.diag(bunch_cov_matrix)), 2))
print(np.round(1000.0 * np.sqrt(np.diag(envelope.cov_matrix)), 2))
