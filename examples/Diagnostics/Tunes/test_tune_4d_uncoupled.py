"""Test one-turn tune estimation in uncoupled lattice using 4D normalization.

See comments on `test_tune.py`. This example is the same, but the tunes are
estimated using the `BunchTuneAnalysis4D` class.  Since there is no coupling
in the lattice, the (eigen)tunes {nu1, nu2} should be the same as horizontal
and vertical tunes {nux, nuy}.
"""
import math
import os
import pathlib
import random
from pprint import pprint

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTuneAnalysis
from orbit.core.bunch import BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist2D
from orbit.bunch_generators import WaterBagDist2D
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.utils.consts import mass_proton

from utils import make_lattice


# Setup
# ------------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


# Initialize lattice and bunch
# ------------------------------------------------------------------------------------

lattice = make_lattice()

bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.000)


# Analyze transfer matrix
# ------------------------------------------------------------------------------------

def build_norm_matrix_from_twiss_2d(alpha: float, beta: float) -> np.ndarray:
    norm_matrix_inv = np.array([[beta, 0.0], [-alpha, 1.0]]) * np.sqrt(1.0 / beta)
    norm_matrix = np.linalg.inv(norm_matrix_inv)
    return norm_matrix


matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
lattice_params = matrix_lattice.getRingParametersDict()
pprint(lattice_params)

lattice_alpha_x = lattice_params["alpha x"]
lattice_alpha_y = lattice_params["alpha y"]
lattice_beta_x = lattice_params["beta x [m]"]
lattice_beta_y = lattice_params["beta y [m]"]

norm_matrix = np.zeros((4, 4))
norm_matrix[0:2, 0:2] = build_norm_matrix_from_twiss_2d(lattice_alpha_x, lattice_beta_x)
norm_matrix[2:4, 2:4] = build_norm_matrix_from_twiss_2d(lattice_alpha_y, lattice_beta_y)

print("Normalization matrix V^{-1}:")
print(norm_matrix)


# Add tune diagnostic node
# ------------------------------------------------------------------------------------
    
tune_node = TeapotTuneAnalysisNode()
tune_node.setNormMatrix(norm_matrix)
lattice.getNodes()[0].addChildNode(tune_node, 0)


# Generate phase space distribution
# ------------------------------------------------------------------------------------

bunch_emitt_x = 25.0e-06
bunch_emitt_y = 25.0e-06
bunch_twiss_x = TwissContainer(lattice_alpha_x, lattice_beta_x, bunch_emitt_x)
bunch_twiss_y = TwissContainer(lattice_alpha_y, lattice_beta_y, bunch_emitt_y)
bunch_dist_xy = GaussDist2D(bunch_twiss_x, bunch_twiss_y)

n_parts = 1000
for index in range(n_parts):
    (x, xp, y, yp) = bunch_dist_xy.getCoordinates()
    z = random.uniform(-25.0, 25.0)
    dE = 0.0
    bunch.addParticle(x, xp, y, yp, z, dE)


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

filename = "bunch.dat"
filename = os.path.join(output_dir, filename)
bunch.dumpBunch(filename)


# Analysis    
# ------------------------------------------------------------------------------------

# Collect phase data from bunch
phase_data = {}
for i in range(bunch.getSize()):
    data = tune_node.getData(bunch, i)
    for key in data:
        if key in phase_data:
            phase_data[key].append(data[key])
        else:
            phase_data[key] = []

phase_data = pd.DataFrame(phase_data)
print(phase_data)

# Read phase data from file
particles = np.loadtxt(filename, comments="%")
particles = pd.DataFrame(
    particles, 
    columns=[  # https://github.com/PyORBIT-Collaboration/PyORBIT3/issues/78
        "x", 
        "xp", 
        "y", 
        "yp", 
        "z",
        "dE",   
        "phase_1",
        "phase_2",
        "tune_1",
        "tune_2",
        "action_1",
        "action_2",
    ] 
)
print(particles.iloc[:, 6:])

# Check against tune from transfer matrix
tune_x_true = lattice_params["fractional tune x"]
tune_y_true = lattice_params["fractional tune y"]
tune_x_pred = np.mean(phase_data["tune_x"])
tune_y_pred = np.mean(phase_data["tune_y"])

tune_x_err = tune_x_pred - tune_x_true
tune_y_err = tune_y_pred - tune_y_true

print("tune_x_err", tune_x_err)
print("tune_y_err", tune_y_err)

assert np.abs(tune_x_err) < 1.00e-08
assert np.abs(tune_y_err) < 1.00e-08