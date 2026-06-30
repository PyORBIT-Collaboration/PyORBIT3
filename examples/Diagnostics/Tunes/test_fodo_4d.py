"""Test one-turn tune estimation in coupled lattice."""

import argparse
import math
import os
import pathlib

import numpy as np
import pandas as pd

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.diagnostics.matrix import TransferMatrixAnalysis
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.teapot import SolenoidTEAPOT
from orbit.utils.consts import mass_proton

from utils import make_lattice_fodo
from utils import get_tmat


# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--coupled", type=int, default=1)
parser.add_argument("--norm-from", type=str, default="tmat", choices=["tmat", "cov"])
args = parser.parse_args()

# Setup
path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)

# Initialize lattice and bunch
lattice = make_lattice_fodo()
if args.coupled:
    node = SolenoidTEAPOT(length=0.5, B=0.25)
    node.setUsageFringeFieldIN(False)
    node.setUsageFringeFieldOUT(False)
    lattice.addNode(node)
lattice.initialize()

bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.000)

# Calculate transfer matrix
M = get_tmat(lattice, bunch)
tmat_analysis = TransferMatrixAnalysis(M)
tune_1_true = tmat_analysis.eigtunes[0]
tune_2_true = tmat_analysis.eigtunes[1]

# Calculate matched covariance matrix
eps_1 = 0.1e-06  # mode 1 rms emittance
eps_2 = 0.1e-06  # mode 2 rms emittance
S_matched = tmat_analysis.cov_matrix(eps_1, eps_2)

# Add tune diagnostic node
tune_node = TeapotTuneAnalysisNode()
if args.norm_from == "tmat":
    tune_node.setNormMatrixFromTransferMatrix(M)
else:
    tune_node.setNormMatrixFromCovMatrix(S_matched)
lattice.getNodes()[0].addChildNode(tune_node, 0)

# Generate particles

rng = np.random.default_rng(123)
particles = np.zeros((1000, 6))
particles[:, :4] = rng.multivariate_normal(
    mean=np.zeros(4), cov=S_matched, size=particles.shape[0]
)
particles[:, 4] = rng.uniform(-25.0, 25.0, size=particles.shape[0])
particles[:, 5] = 0.0

for index in range(particles.shape[0]):
    bunch.addParticle(*particles[index])

# Track particles
for turn in range(10):
    lattice.trackBunch(bunch)

    twiss_calc = BunchTwissAnalysis()
    twiss_calc.analyzeBunch(bunch)
    xrms = 1000.0 * math.sqrt(twiss_calc.getCorrelation(0, 0))
    yrms = 1000.0 * math.sqrt(twiss_calc.getCorrelation(2, 2))
    print("turn={} xrms={:0.3f} yrms={:0.3f}".format(turn + 1, xrms, yrms))

# Analysis
phase_data = tune_node.getData(bunch)  # phase_data = pd.DataFrame(phase_data)
print(phase_data)

tune_1_calc = np.mean(phase_data["tune_1"])
tune_2_calc = np.mean(phase_data["tune_2"])
tune_1_err = tune_1_calc - tune_1_true
tune_2_err = tune_2_calc - tune_2_true

print("tune_1_true", tune_1_true)
print("tune_1_calc", tune_1_calc)
print("tune_2_true", tune_2_true)
print("tune_2_calc", tune_2_calc)
print("tune_1_err", tune_1_err)
print("tune_2_err", tune_2_err)

assert np.abs(tune_1_err) < 1.00e-08
assert np.abs(tune_2_err) < 1.00e-08
