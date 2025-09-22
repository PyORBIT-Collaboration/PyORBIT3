import argparse
import os
import pathlib

import numpy as np
import matplotlib.pyplot as plt

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch
from orbit.diagnostics import BunchHistogram
from orbit.diagnostics import BunchHistogram1D
from orbit.diagnostics import BunchHistogram2D
from orbit.diagnostics import BunchHistogram3D
from orbit.bunch_utils import collect_bunch

plt.style.use("style.mplstyle")


# Parse args
parser = argparse.ArgumentParser()
parser.add_argument("--nsamp", type=int, default=10_000)
parser.add_argument("--nbins", type=int, default=64)
args = parser.parse_args()

# Make output directory
path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)

# Setup MPI
_mpi_comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
_mpi_rank = orbit_mpi.MPI_Comm_rank(_mpi_comm)

# Generate particle distribution
rng = np.random.default_rng(123)
X = rng.normal(size=(args.nsamp, 2))
X = X / np.linalg.norm(X, axis=1)[:, None]
X[:, 0] -= 0.75 * X[:, 1] ** 2
X = X + rng.normal(size=X.shape, scale=0.25)
X = X / np.std(X, axis=0)
X = np.hstack([X, np.zeros((X.shape[0], 4))])

# Create bunch
bunch = Bunch()
for i in range(X.shape[0]):
    bunch.addParticle(*X[i, :])

# Create histogram diagnostic
grid_limits = [(-4.0, 4.0)]
grid_shape = (args.nbins,)

axis = (0,)
diag = BunchHistogram(
    axis=axis,
    shape=grid_shape,
    limits=grid_limits,
    method=None,
    normalize=True,
    output_dir=output_dir,
)

# Compute histogram
diag.track(bunch)

# Plot histogram
if _mpi_rank == 0:
    fig, ax = plt.subplots(figsize=(5, 2))
    ax.plot(diag.grid_coords[0], diag.grid_values, lw=1.5, color=None, label="PyORBIT")
    grid_values_np, _ = np.histogram(X[:, axis], bins=diag.grid_edges[0], density=True)
    ax.plot(diag.grid_coords[0], grid_values_np, lw=1.5, color=None, label="NumPy")
    ax.legend()
    
    filename = os.path.join(output_dir, "fig_hist")
    plt.savefig(filename)
    plt.show()
