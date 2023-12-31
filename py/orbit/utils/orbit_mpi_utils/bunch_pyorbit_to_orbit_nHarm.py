import math

from orbit.core.bunch import Bunch

# pyORBIT MPI module import
from orbit.core import orbit_mpi
from orbit.core.orbit_mpi import mpi_datatype
from orbit.core.orbit_mpi import mpi_op


def bunch_pyorbit_to_orbit_nHarm(ringLength, nHarm, pyOrbitBunch, name_of_orbit_mpi_bunch_file):
    """
    Translates pyORBIT bunch to ORBIT_MPI bunch, incorporating RF
    harmonic number, and dumps it into a file.
    The ring length should be defined in the input (in meters).
    Lines in bunch files:
    ORBIT_MPI: x[mm] xp[mrad] y[mm] yp[mrad] phi[rad] dE[GeV].
    pyORBIT:   x[m]  xp[rad]  y[m]  yp[rad]  z[m]     dE[GeV]
    """
    pi2 = 2.0 * math.pi
    zfac = pi2 * nHarm / ringLength
    b = pyOrbitBunch
    # Take the MPI Communicator from bunch: it could be different
    # from MPI_COMM_WORLD
    comm = pyOrbitBunch.getMPIComm()
    rank = orbit_mpi.MPI_Comm_rank(comm)
    size = orbit_mpi.MPI_Comm_size(comm)
    main_rank = 0

    # n_parts_arr - array of size number of CPUs,
    # Contains the number of macroparticles on each CPU
    n_parts_arr = [0] * size
    n_parts_arr[rank] = b.getSize()
    n_parts_arr = orbit_mpi.MPI_Allreduce(n_parts_arr, mpi_datatype.MPI_INT, mpi_op.MPI_SUM, comm)

    file_out = None
    if rank == main_rank:
        file_out = open(name_of_orbit_mpi_bunch_file, "w")

    if rank == main_rank:
        for i in range(n_parts_arr[rank]):
            x = b.x(i) * 1000.0
            px = b.px(i) * 1000.0
            y = b.y(i) * 1000.0
            py = b.py(i) * 1000.0
            z = -(math.fmod(b.z(i) * zfac, pi2))
            if z > math.pi:
                z = z - 2 * math.pi
            if z < -math.pi:
                z = z + 2 * math.pi
            dE = b.dE(i)
            file_out.write(str(x) + " " + str(px) + " " + str(y) + " " + str(py) + " " + str(z) + " " + str(dE) + "\n")

    # MPI_Barrier command is not necessary.
    orbit_mpi.MPI_Barrier(comm)

    val_arr = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    for i_cpu in range(1, size):
        # Again, MPI_Barrier command is not necessary.
        orbit_mpi.MPI_Barrier(comm)
        for i in range(n_parts_arr[i_cpu]):
            if rank == main_rank:
                # Get the coordinate array
                (x, px, y, py, z, dE) = orbit_mpi.MPI_Recv(mpi_datatype.MPI_DOUBLE, i_cpu, 222, comm)
                file_out.write(str(x) + " " + str(px) + " " + str(y) + " " + str(py) + " " + str(z) + " " + str(dE) + "\n")
            elif rank == i_cpu:
                # send the coordinate array
                x = b.x(i) * 1000.0
                px = b.px(i) * 1000.0
                y = b.y(i) * 1000.0
                py = b.py(i) * 1000.0
                z = -(math.fmod(b.z(i) * zfac, pi2))
                if z > math.pi:
                    z = z - 2 * math.pi
                if z < -math.pi:
                    z = z + 2 * math.pi
                dE = b.dE(i)
                val_arr = (x, px, y, py, z, dE)
                orbit_mpi.MPI_Send(val_arr, mpi_datatype.MPI_DOUBLE, main_rank, 222, comm)

    if rank == main_rank:
        file_out.close()
