"""Test one-turn tune estimation in uncoupled lattice.

This example tracks a Gaussian distribution through a FODO lattice. The tunes
are estimated from the phase space coordinates before/after tracking using the
`BunchTuneAnalysis` class.
"""
import argparse
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
from orbit.diagnostics.matrix import TransferMatrixAnalysis
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_lattice_sns
from utils import get_tmat


# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--norm-from", type=str, default="tmat", choices=["tmat", "cov", "twiss"])
args = parser.parse_args()

# Setup
path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)

# Initialize lattice and bunch
lattice = make_lattice_sns()

bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.000)

# Calculate transfer matrix
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

# Add tune diagnostic node
tune_node = TeapotTuneAnalysisNode()

if args.norm_from == "twiss":
    tune_node.setNormMatrixFromTwiss(
        betax=lattice_beta_x,
        alphax=lattice_alpha_x,
        etax=lattice_eta_x,
        etapx=lattice_etap_x,
        betay=lattice_beta_y,
        alphay=lattice_alpha_y,
    )
elif args.norm_from == "tmat":
    tmat = get_tmat(lattice, bunch)
    tune_node.setNormMatrixFromTransferMatrix(tmat)
elif args.norm_from == "cov":
    tmat = get_tmat(lattice, bunch)
    tmat_analysis = TransferMatrixAnalysis(tmat)
    cov_matrix = tmat_analysis.cov_matrix(1e-7, 1e-7)
    tune_node.setNormMatrixFromCovMatrix(cov_matrix)
else:
    raise ValueError("Invalid norm_from argument")

lattice.getNodes()[0].addChildNode(tune_node, 0)

# Generate particles
emittance_x = 0.1e-06
emittance_y = 0.1e-06
twiss_x = TwissContainer(lattice_alpha_x, lattice_beta_x, emittance_x)
twiss_y = TwissContainer(lattice_alpha_y, lattice_beta_y, emittance_y)
dist = GaussDist2D(twiss_x, twiss_y)

for index in range(1000):
    x, xp, y, yp = dist.getCoordinates()
    z = random.uniform(-25.0, 25.0)
    bunch.addParticle(x, xp, y, yp, z, 0.0)

# Track particles
for turn in range(10):
    lattice.trackBunch(bunch)

    twiss_calc = BunchTwissAnalysis()
    twiss_calc.analyzeBunch(bunch)
    xrms = 1000.0 * math.sqrt(twiss_calc.getCorrelation(0, 0))
    yrms = 1000.0 * math.sqrt(twiss_calc.getCorrelation(2, 2))
    print("turn={} xrms={:0.3f} yrms={:0.3f}".format(turn + 1, xrms, yrms))

# Test writing to file
filename = "bunch.dat"
filename = os.path.join(output_dir, filename)
bunch.dumpBunch(filename)

# Collect phase data from bunch
phase_data = tune_node.getData(bunch)
phase_data = pd.DataFrame(phase_data)

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
    ],
)
print(particles.iloc[:, 6:])

# Check against tune from transfer matrix
tune_x_true = lattice_params["fractional tune x"]
tune_y_true = lattice_params["fractional tune y"]
tune_x_calc = np.mean(phase_data["tune_1"])
tune_y_calc = np.mean(phase_data["tune_2"])

tune_x_err = tune_x_calc - tune_x_true
tune_y_err = tune_y_calc - tune_y_true

print("tune_x_true", tune_x_true)
print("tune_x_calc", tune_x_calc)
print("tune_y_true", tune_y_true)
print("tune_y_calc", tune_y_calc)
print("tune_x_err", tune_x_err)
print("tune_y_err", tune_y_err)

assert np.abs(tune_x_err) < 1.00e-06
assert np.abs(tune_y_err) < 1.00e-06
