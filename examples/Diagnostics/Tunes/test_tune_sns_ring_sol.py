import math
import os
import pathlib
import random
from pprint import pprint

import numpy as np
import pandas as pd

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist2D
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.teapot import TEAPOT_Ring
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.teapot import teapot
from orbit.utils.consts import mass_proton

# Setup
# ------------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)

# Initialize lattice and bunch
# ------------------------------------------------------------------------------------

lattice = TEAPOT_Ring()
lattice.readMADX("inputs/sns_ring.lat", "rnginjsol")
lattice.initialize()

for node in lattice.getNodes():
    # Turn off fringe fields.
    try:
        node.setUsageFringeFieldIN(False)
        node.setUsageFringeFieldOUT(False)
    except:
        pass

    # Turn off solenoid.
    for name in ["scbdsol_c13a", "scbdsol_c13b"]:
        node = lattice.getNodeForName(name)
        node.setParam("B", 0.15 * 0.5)

bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.000)

# Analyze transfer matrix
# ------------------------------------------------------------------------------------

from orbit.diagnostics.eig import calc_eigtune, normalize_eigvec, calc_norm_matrix_from_eigvecs

# Estimate transfer matrix
matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)

M = np.zeros((4, 4))
for i in range(4):
    for j in range(4):
        M[i, j] = matrix_lattice.getOneTurnMatrix().get(i, j)

# Calculate eigenvalues and eigenvectors
eigvals, eigvecs = np.linalg.eig(M)
eigvals = eigvals[::2]
eigvecs = eigvecs[:, ::2].T

for i in range(eigvecs.shape[0]):
    eigvecs[i] = normalize_eigvec(eigvecs[i])

# Calculate tunes from transfer matrix
tune_1_true = calc_eigtune(eigvals[0])
tune_2_true = calc_eigtune(eigvals[1])

# Calculate normalization matrix from transfer matrix
V_inv = calc_norm_matrix_from_eigvecs(*eigvecs)
V = np.linalg.inv(V_inv)

# Print normalization matrix
print("Normalization matrix V^{-1}:")
print(V_inv)

# Add tune diagnostic node
# ------------------------------------------------------------------------------------

tune_node = TeapotTuneAnalysisNode()
tune_node.setNormMatrix(V_inv)
lattice.getNodes()[0].addChildNode(tune_node, 0)

# Generate phase space distribution
# ------------------------------------------------------------------------------------

rng = np.random.default_rng()

n = 1000
eps_1 = 0.1e-06  # mode 1 rms emittance
eps_2 = 0.1e-06  # mode 2 rms emittance

# Generate particles in normalized phase space
particles = np.zeros((n, 6))
particles[:, (0, 1)] = rng.normal(size=(n, 2), scale=np.sqrt(eps_1))
particles[:, (2, 3)] = rng.normal(size=(n, 2), scale=np.sqrt(eps_2))
particles[:, 4] = rng.uniform(-25.0, 25.0, size=n)
particles[:, 5] = 0.0

# Unnormalize transverse coordinates (match to lattice)
particles[:, :4] = np.matmul(particles[:, :4], V.T)

# Add particles to bunch
for index in range(n):
    bunch.addParticle(*particles[index])

# Tracking
# ------------------------------------------------------------------------------------

n_turns = 10
for turn in range(n_turns):
    lattice.trackBunch(bunch)

    twiss_calc = BunchTwissAnalysis()
    twiss_calc.analyzeBunch(bunch)
    xrms = math.sqrt(twiss_calc.getCorrelation(0, 0)) * 1000.0
    yrms = math.sqrt(twiss_calc.getCorrelation(2, 2)) * 1000.0
    print("turn={} xrms={:0.3f} yrms={:0.3f}".format(turn + 1, xrms, yrms))

# Analysis
# ------------------------------------------------------------------------------------

# Collect phase data from bunch
phase_data = tune_node.getData(bunch)
phase_data = pd.DataFrame(phase_data)
print(phase_data)

# Check average tune vs. transfer matrix
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

assert np.abs(tune_1_err) < 1.00e-06
assert np.abs(tune_2_err) < 1.00e-06
