from orbit.core.bunch import Bunch
from orbit.core import orbit_mpi

import numpy as np
from numpy.typing import NDArray


def collect_bunch(
    bunch: Bunch,
) -> dict[str, np.float64 | np.int32 | NDArray[np.float64]]:
    """Collects attributes from a PyOrbit Bunch across all MPI ranks and returns it as a dictionary.
    Parameters
    ----------
    bunch : Bunch
        The PyOrbit::Bunch object from which to collect attributes.
    Returns
    -------
    dict[str, np.float64 | np.int32 | NDArray[np.float64]]:
        By default this returns a dictionary containing the following keys and their corresponding values:
        - "x": particle x-coordinates [m]
        - "xp": particle x-momenta [rad]
        - "y": particle y-coordinates [m]
        - "yp": particle y-momenta [rad]
        - "z": particle longitudinal coordinates [m]
        - "dE": particle energy deviations [GeV]
        - "sync_part_coords": coordinates of the synchronous particle (x,y,z) [m]
        - "sync_part_kin_energy": kinetic energy of the synchronous particle [GeV]
        - "sync_part_momentum": momentum of the synchronous particle  [GeV/c]
        - "sync_part_beta": beta of the synchronous particle
        - "sync_part_gamma": gamma of the synchronous particle
        - "sync_part_time": time of the synchronous particle [s]
        - Any additional attributes defined in the bunch.
    """
    n_particles = bunch.getSize()

    if n_particles == 0:
        return {}

    mpi_comm = bunch.getMPIComm()  # orbit_mpi.mpi_comm.MPI_COMM_WORLD
    mpi_rank = orbit_mpi.MPI_Comm_rank(mpi_comm)
    mpi_size = orbit_mpi.MPI_Comm_size(mpi_comm)

    if mpi_rank == 0:
        bunch_dict = {"x": [], "xp": [], "y": [], "yp": [], "z": [], "dE": []}
        for attr in bunch.bunchAttrDoubleNames():
            bunch_dict[attr] = np.float64(bunch.bunchAttrDouble(attr))

        for attr in bunch.bunchAttrIntNames():
            bunch_dict[attr] = np.int32(bunch.bunchAttrInt(attr))

        sync_part = bunch.getSyncParticle()

        bunch_dict |= {
            "sync_part_coords": np.array(sync_part.pVector()),
            "sync_part_kin_energy": np.float64(sync_part.kinEnergy()),
            "sync_part_momentum": np.float64(sync_part.momentum()),
            "sync_part_beta": np.float64(sync_part.beta()),
            "sync_part_gamma": np.float64(sync_part.gamma()),
            "sync_part_time": np.float64(sync_part.time()),
        }

        for i in range(n_particles):
            bunch_dict["x"].append(bunch.x(i))
            bunch_dict["xp"].append(bunch.xp(i))
            bunch_dict["y"].append(bunch.y(i))
            bunch_dict["yp"].append(bunch.yp(i))
            bunch_dict["z"].append(bunch.z(i))
            bunch_dict["dE"].append(bunch.dE(i))

    mpi_tag = 42  # not sure this is necessary; seems like it can be any integer
    for cpu_idx in range(1, mpi_size):
        for i in range(n_particles):
            if mpi_rank == cpu_idx:
                coord_arr = (
                    bunch.x(i),
                    bunch.xp(i),
                    bunch.y(i),
                    bunch.yp(i),
                    bunch.z(i),
                    bunch.dE(i),
                )
                orbit_mpi.MPI_Send(
                    coord_arr, orbit_mpi.mpi_datatype.MPI_DOUBLE, 0, mpi_tag, mpi_comm
                )
            elif mpi_rank == 0:
                coord_arr = orbit_mpi.MPI_Recv(
                    orbit_mpi.mpi_datatype.MPI_DOUBLE, cpu_idx, mpi_tag, mpi_comm
                )
                bunch_dict["x"].append(coord_arr[0])
                bunch_dict["xp"].append(coord_arr[1])
                bunch_dict["y"].append(coord_arr[2])
                bunch_dict["yp"].append(coord_arr[3])
                bunch_dict["z"].append(coord_arr[4])
                bunch_dict["dE"].append(coord_arr[5])

    if mpi_rank == 0:
        for k, v in bunch_dict.items():
            if isinstance(v, list):
                bunch_dict[k] = np.array(v, dtype=np.float64)
        return bunch_dict
