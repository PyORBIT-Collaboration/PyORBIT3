"""Benchmark Danilov envelope tracker vs. PIC for test particles."""

import argparse
import copy
import os
import pathlib
from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from tqdm import tqdm
from tqdm import trange

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.envelope import DanilovEnvelope
from orbit.envelope import DanilovEnvelopeTracker
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.space_charge.sc2p5d import SC2p5D_AccNode
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_fodo_lattice
from utils import BunchMonitor
from utils import rms_ellipse_params

plt.style.use("style.mplstyle")


# Parse arguments
# --------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--phase-adv", type=float, default=100.0)
parser.add_argument("--intensity", type=float, default=7.0e14)
parser.add_argument("--eps_x", type=float, default=10.00e-06)
parser.add_argument("--eps_y", type=float, default=10.00e-06)
parser.add_argument("--max-part-length", type=float, default=0.1)
parser.add_argument("--mismatch", type=float, default=0.0)
parser.add_argument("--periods", type=int, default=100)
args = parser.parse_args()


# Setup
# --------------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


# Set up simulation
# --------------------------------------------------------------------------------------

envelope = DanilovEnvelope(
    eps_1=(args.eps_x + args.eps_y),
    eps_2=0.0,
    mass=0.938,
    kin_energy=1.0,
    length=100.0,
    line_density=(args.intensity / 100.0),
    params=None,
)

lattice = make_fodo_lattice(
    phase_adv_x=np.radians(args.phase_adv),
    phase_adv_y=np.radians(args.phase_adv),
    length=5.0,
    mass=envelope.mass,
    kin_energy=envelope.kin_energy,
    max_part_length=args.max_part_length,
    verbose=1,
)

tracker = DanilovEnvelopeTracker(lattice, path_length_max=args.max_part_length)
tracker.match_zero_sc(envelope, method="2d")
tracker.match(envelope, method="replace_avg", periods_avg=20, verbose=2)
envelope_init = envelope.copy()


# Set test particle coordinates
# --------------------------------------------------------------------------------------

test_particles = [
    [0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.005, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.010, 0.0, 0.0, 0.0, 0.0, 0.0],
]
test_particles = np.array(test_particles)


# Track envelope with test particles
# --------------------------------------------------------------------------------------
    
particles_tbt = {}
particles_tbt["envelope"] = []
particles_tbt["bunch"] = []

envelope = envelope_init.copy()
particles = test_particles.copy()
for period in range(args.periods):
    envelope, particles = tracker.track_particles(envelope, particles=particles)

    cov_matrix = envelope.cov()
    cov_matrix *= 1.00e+06
    xrms = np.sqrt(cov_matrix[0, 0])
    yrms = np.sqrt(cov_matrix[2, 2])
    print(f"turn={period} xrms={xrms:0.3f} yrms={yrms:0.3f}")
    
    particles_tbt["envelope"].append(particles.copy())
    

# Track bunch with test particles
# --------------------------------------------------------------------------------------

def set_particle_macrosizes(bunch: Bunch, macrosizes: list[float]) -> Bunch:
    bunch.addPartAttr("macrosize")  # sets macrosize=0 for all particles
    attribute_array_index = 0
    for index in range(bunch.getSize()):
        bunch.partAttrValue("macrosize", index, attribute_array_index, macrosizes[index])
    return bunch


lattice = make_fodo_lattice(
    phase_adv_x=np.radians(args.phase_adv),
    phase_adv_y=np.radians(args.phase_adv),
    length=5.0,
    mass=envelope.mass,
    kin_energy=envelope.kin_energy,
    max_part_length=args.max_part_length,
    verbose=1,
)

sc_calc = SpaceChargeCalc2p5D(64, 64, 1)
sc_path_length_min = 1.00e-06
sc_nodes = setSC2p5DAccNodes(lattice, sc_path_length_min, sc_calc)

bunch_size = 128_000
bunch = envelope_init.to_bunch(env=False, size=bunch_size)

for i in range(test_particles.shape[0]):
    bunch.addParticle(*test_particles[i])

macrosize = args.intensity / bunch_size
macrosizes = macrosize * np.ones(bunch.getSize())
macrosizes[bunch_size:] = 0
set_particle_macrosizes(bunch, macrosizes)

particles_tbt["bunch"] = []
for period in range(args.periods):
    lattice.trackBunch(bunch)

    bunch_calc = BunchTwissAnalysis()
    bunch_calc.analyzeBunch(bunch)
    xrms = 1000.0 * np.sqrt(bunch_calc.getCorrelation(0, 0))
    yrms = 1000.0 * np.sqrt(bunch_calc.getCorrelation(2, 2))
    print(f"turn={period} xrms={xrms:0.3f} yrms={yrms:0.3f}")

    particles = []
    for i in range(bunch.getSize() - test_particles.shape[0], bunch.getSize()):
        particles.append([bunch.x(i), bunch.xp(i), bunch.y(i), bunch.yp(i), bunch.z(i), bunch.dE(i)])
    particles = np.array(particles)
    particles_tbt["bunch"].append(particles.copy())


# Analysis
# --------------------------------------------------------------------------------------

# Collect data
for key in ["envelope", "bunch"]:
    particles_tbt[key] = np.stack(particles_tbt[key])

# Plot x-x'
fig, axs = plt.subplots(ncols=2, figsize=(6.0, 3.0), sharex=True, sharey=True)
for ax, key in zip(axs, ["envelope", "bunch"]):
    particles = particles_tbt[key] * 1000.0
    for index in range(particles.shape[1]):
        ax.scatter(particles[:, index, 0], particles[:, index, 1], s=3)

cov_matrix = envelope.cov()
cov_matrix = cov_matrix[0:2, 0:2]
cov_matrix = cov_matrix * 1.00e+06
cx, cy, angle = rms_ellipse_params(cov_matrix)
for ax in axs:
    ax.add_patch(
        patches.Ellipse(
            xy=(0.0, 0.0),
            width=(4.0 * cx),
            height=(4.0 * cy),
            angle=(-np.degrees(angle)),
            color="black",
            ec="none",
            alpha=0.1,
            zorder=0,
        )
    )
for ax in axs:
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("xp [mrad]")
    
filename = "fig_benchmark_particles_xxp.png"
filename = os.path.join(output_dir, filename)
plt.savefig(filename, dpi=300)
plt.show()


