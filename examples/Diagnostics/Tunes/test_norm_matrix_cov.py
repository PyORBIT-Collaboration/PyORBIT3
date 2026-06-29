import math
import os
import pathlib

import numpy as np
import pandas as pd

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.teapot import TEAPOT_Ring
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from orbit.diagnostics.eig import calc_eigtune
from orbit.diagnostics.eig import normalize_eigvec
from orbit.diagnostics.eig import calc_norm_matrix_from_eigvecs
from orbit.diagnostics.eig import calc_norm_matrix_from_cov
from orbit.diagnostics.eig import calc_norm_matrix_from_tmat

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
    try:
        node.setUsageFringeFieldIN(False)
        node.setUsageFringeFieldOUT(False)
    except:
        pass

    for name in ["scbdsol_c13a", "scbdsol_c13b"]:
        node = lattice.getNodeForName(name)
        node.setParam("B", 0.15 * 0.5)

bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.000)

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)

M = np.zeros((4, 4))
for i in range(4):
    for j in range(4):
        M[i, j] = matrix_lattice.getOneTurnMatrix().get(i, j)

tune_node = TeapotTuneAnalysisNode()
lattice.getNodes()[0].addChildNode(tune_node, 0)

rng = np.random.default_rng()

V_inv = calc_norm_matrix_from_tmat(M)
V = np.linalg.inv(V_inv)

# From transfer matrix
# ------------------------------------------------------------------------------------

tune_node.setNormMatrixFromTransferMatrix(M)

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

n_turns = 10
for turn in range(n_turns):
    lattice.trackBunch(bunch)

    twiss_calc = BunchTwissAnalysis()
    twiss_calc.analyzeBunch(bunch)
    xrms = math.sqrt(twiss_calc.getCorrelation(0, 0)) * 1000.0
    yrms = math.sqrt(twiss_calc.getCorrelation(2, 2)) * 1000.0

# Collect phase data from bunch
phase_data = tune_node.getData(bunch)
phase_data = pd.DataFrame(phase_data)

# Check average tune vs. transfer matrix
tune_1_calc = np.mean(phase_data["tune_1"])
tune_2_calc = np.mean(phase_data["tune_2"])

print("tune_1_calc", tune_1_calc)
print("tune_2_calc", tune_2_calc)


# From transfer matrix
# ------------------------------------------------------------------------------------

S = np.matmul(V, V.T)
tune_node.setNormMatrixFromCovMatrix(S)

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
bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.000)
for index in range(n):
    bunch.addParticle(*particles[index])

n_turns = 10
for turn in range(n_turns):
    lattice.trackBunch(bunch)

    twiss_calc = BunchTwissAnalysis()
    twiss_calc.analyzeBunch(bunch)
    xrms = math.sqrt(twiss_calc.getCorrelation(0, 0)) * 1000.0
    yrms = math.sqrt(twiss_calc.getCorrelation(2, 2)) * 1000.0

# Collect phase data from bunch
phase_data = tune_node.getData(bunch)
phase_data = pd.DataFrame(phase_data)

# Check average tune vs. transfer matrix
tune_1_calc = np.mean(phase_data["tune_1"])
tune_2_calc = np.mean(phase_data["tune_2"])

print("tune_1_calc", tune_1_calc)
print("tune_2_calc", tune_2_calc)