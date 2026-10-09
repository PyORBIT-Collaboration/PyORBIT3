"""Test updating normalization matrix during tracking.

Average tunes are printed out after each turn. The normalization matrix
is updated twice during tracking. A warning message should print that
indicates the normalization matrix has changed, and that the tunes
will be inaccurate until the next tracked turn.
"""
import os
import pathlib
import random
from pprint import pprint

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist2D
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_lattice_sns

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

tune_node.setNormMatrixFromTwiss(
    betax=lattice_beta_x,
    alphax=lattice_alpha_x,
    etax=lattice_eta_x,
    etapx=lattice_etap_x,
    betay=lattice_beta_y,
    alphay=lattice_alpha_y,
)

lattice.getNodes()[0].addChildNode(tune_node, 0)

# Generate particles
emittance_x = 0.1e-06
emittance_y = 0.1e-06
twiss_x = TwissContainer(lattice_alpha_x, lattice_beta_x, emittance_x)
twiss_y = TwissContainer(lattice_alpha_y, lattice_beta_y, emittance_y)
dist = GaussDist2D(twiss_x, twiss_y)

for index in range(10_000):
    x, xp, y, yp = dist.getCoordinates()
    z = random.uniform(-25.0, 25.0)
    bunch.addParticle(x, xp, y, yp, z, 0.0)

# Track particles
for turn in range(20):
    lattice.trackBunch(bunch)


    if turn >= 0:
        phase_data = tune_node.getData(bunch)
        nu_x = np.mean(phase_data["tune_1"])
        nu_y = np.mean(phase_data["tune_2"])

        print("turn{} nux={:0.5f} nuy={:0.5f}".format(turn, nu_x, nu_y))

    if turn == 5:
        tune_node.setNormMatrixFromBunch(bunch)

    if turn == 10:
        tune_node.setNormMatrixFromTwiss(
            betax=lattice_beta_x,
            alphax=lattice_alpha_x,
            etax=lattice_eta_x,
            etapx=lattice_etap_x,
            betay=lattice_beta_y,
            alphay=lattice_alpha_y,
        )
