"""Test one-turn tune estimation in uncoupled lattice.

This example tracks a Gaussian distribution through a FODO lattice. The tunes
are estimated from the phase space coordinates before/after tracking using the
`BunchTuneAnalysis` class.
"""
import math
import os
import pathlib
import random
from pprint import pprint

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist2D
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_lattice


# Setup
# ------------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


# Lattice
# ------------------------------------------------------------------------------------

lattice = make_lattice()

bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.000)

# Compute lattice parameters from one-turn transfer matrix
matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
lattice_params = matrix_lattice.getRingParametersDict()
pprint(lattice_params)

# Store some parameters as variables
lattice_alpha_x = lattice_params["alpha x"]
lattice_alpha_y = lattice_params["alpha y"]
lattice_beta_x = lattice_params["beta x [m]"]
lattice_beta_y = lattice_params["beta y [m]"]
lattice_eta_x = lattice_params["dispersion x [m]"]
lattice_etap_x = lattice_params["dispersion prime x"]


# Tune diagnostics node
# ------------------------------------------------------------------------------------

tune_node = TeapotTuneAnalysisNode()
tune_node.assignTwiss(
    betax=lattice_beta_x, 
    alphax=lattice_alpha_x, 
    etax=lattice_eta_x, 
    etapx=lattice_etap_x, 
    betay=lattice_beta_y,
    alphay=lattice_alpha_y,
)
lattice.getNodes()[0].addChildNode(tune_node, 0)


# Bunch
# ------------------------------------------------------------------------------------

# Generate a matched transverse phase space distribution. The longitudinal 
# distribution will be uniform in position (z) and a delta function in energy 
# deviation (dE).
emittance_x = 25.0e-06
emittance_y = 25.0e-06
bunch_twiss_x = TwissContainer(lattice_alpha_x, lattice_beta_x, emittance_x)
bunch_twiss_y = TwissContainer(lattice_alpha_y, lattice_beta_y, emittance_y)
bunch_dist = GaussDist2D(bunch_twiss_x, bunch_twiss_y)

n_parts = 1000
for index in range(n_parts):
    (x, xp, y, yp) = bunch_dist.getCoordinates()
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
        "phase_x",
        "phase_y",
        "tune_x",
        "tune_y",
        "action_x",
        "action_y",
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