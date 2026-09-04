from orbit.core.bunch import Bunch
from orbit.core import orbit_mpi

import numpy as np


def bunch_from_shared_numpy(global_coords: np.ndarray) -> Bunch:
    """Construct a bunch from the local partition of a shared NumPy array.

    The global particle coordinate array is divided into contiguous, approximately
    equal partitions using the rank and size of ``MPI_COMM_WORLD``. Each rank reads
    its partition directly from ``global_coords``; no particle data is
    transferred through MPI. Consequently, every rank must be able to access
    the complete input array, for example through a shared memory-mapped file.

    Parameters
    ----------
    global_coords : numpy.ndarray
        Global particle coordinates with shape ``(n_particles, 6)`` ordered as
        ``(x, xp, y, yp, z, dE)``.

    Returns
    -------
    Bunch
        A bunch containing the calling rank's contiguous partition of the
        global particle array.

    Raises
    ------
    ValueError
        If ``global_coords`` does not have shape ``(n_particles, 6)``.

    Notes
    -----
    If the number of particles is not evenly divisible by the number of MPI
    ranks, the lowest-numbered ranks receive one additional particle.
    """
    comm_world = orbit_mpi.mpi_comm.MPI_COMM_WORLD
    mpi_size = orbit_mpi.MPI_Comm_size(comm_world)
    rank = orbit_mpi.MPI_Comm_rank(comm_world)

    global_size = global_coords.shape[0]
    base = global_size // mpi_size
    remainder = global_size % mpi_size

    local_size = base + (1 if rank < remainder else 0)
    start_row = rank * base + min(rank, remainder)
    stop_row = start_row + local_size

    return Bunch.from_numpy(global_coords[start_row:stop_row, :])
