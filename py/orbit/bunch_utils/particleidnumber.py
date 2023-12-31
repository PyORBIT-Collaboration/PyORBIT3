#!/usr/bin/env python

"""
"""
import math
import random
import sys
from orbit.core.bunch import Bunch
from orbit.core import orbit_mpi
from orbit.core.orbit_mpi import mpi_comm, mpi_datatype, mpi_op


class ParticleIdNumber:
    """
    This routine adds id numbers to particle in a bunch.
    """

    @staticmethod
    def addParticleIdNumbers(b, fixedidnumber=-1, part_ind=0):
        rank = 0
        numprocs = 1

        mpi_init = orbit_mpi.MPI_Initialized()
        comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD

        if mpi_init:
            rank = orbit_mpi.MPI_Comm_rank(comm)
            numprocs = orbit_mpi.MPI_Comm_size(comm)

        nparts_arr_local = []
        for i in range(numprocs):
            nparts_arr_local.append(0)

        nparts_arr_local[rank] = b.getSize()
        data_type = mpi_datatype.MPI_INT
        op = mpi_op.MPI_SUM

        nparts_arr = orbit_mpi.MPI_Allreduce(nparts_arr_local, data_type, op, comm)

        if b.hasPartAttr("ParticleIdNumber") == 0:
            b.addPartAttr("ParticleIdNumber")

        if fixedidnumber >= 0:
            for i in range(part_ind, b.getSize()):
                b.partAttrValue("ParticleIdNumber", i, 0, fixedidnumber)

        else:
            istart = 0
            if rank == 0:
                istart = 0
            else:
                for i in range(rank):
                    istart = istart + nparts_arr[i]

            for i in range(b.getSize()):
                idnumber = istart + i
                b.partAttrValue("ParticleIdNumber", i, 0, idnumber)
