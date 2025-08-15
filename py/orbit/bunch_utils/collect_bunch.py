import os
# import tempfile
from enum import IntEnum
from typing import Optional, TypedDict

from orbit.core.bunch import Bunch
from orbit.core import orbit_mpi

import numpy as np
from numpy.typing import NDArray


class SyncPartDict(TypedDict):
    coords: NDArray[np.float64]
    kin_energy: np.float64
    momentum: np.float64
    beta: np.float64
    gamma: np.float64
    time: np.float64


class BunchDict(TypedDict):
    coords_array: NDArray[np.float64]
    sync_part: SyncPartDict
    attr: dict[str, np.float64 | np.int32]


class BunchCoord(IntEnum):
    X  = 0
    XP = 1
    Y  = 2
    YP = 3
    Z  = 4
    DE = 5


def collect_bunch(
    bunch: Bunch,
    return_memmap: bool = True,
    output_fname: Optional[str] = None,
) -> BunchDict | None:
    """Collects attributes from a PyOrbit Bunch across all MPI ranks and returns it as a dictionary.
    Parameters
    ----------
    bunch : Bunch
        The PyOrbit::Bunch object from which to collect attributes.
    return_memmap : bool, optional
        Return the bunch coordinates as a memory-mapped NumPy array, otherwise the
        entire array is copied into RAM and returned as normal NDArray. Default is True.
    Returns
    -------
    BunchDict | None
        A dictionary containing the collected bunch attributes. Returns None if not on the root MPI rank or if the global bunch size is 0.
        BunchDict structure:
            {
                "coords": NDArray[np.float64] of shape (N, 6) where N is the total number of macroparticles,
                    and the 6 columns correspond to [x, xp, y, yp, z, dE] in units of [m, rad, m, rad, m, eV], respectively.
                "sync_part": {
                    "coords": NDArray[np.float64] of shape (3,),
                    "kin_energy": np.float64,
                    "momentum": np.float64,
                    "beta": np.float64,
                    "gamma": np.float64,
                    "time": np.float64
                },
                "attributes": {
                    <bunch attribute name>: <attribute value (np.float64 or np.int32)>,
                    ...
                }
            }
    """

    global_size = bunch.getSizeGlobal()

    if global_size == 0:
        return None

    mpi_comm = bunch.getMPIComm()
    mpi_rank = orbit_mpi.MPI_Comm_rank(mpi_comm)

    coords_shape = (bunch.getSizeGlobal(), 6)

    local_rows = bunch.getSize()

    # print(f"[DEBUG] Rank {mpi_rank}: start_row={start_row}, stop_row={stop_row}, local_rows={local_rows} bunch.getSize()={bunch.getSize()}")

    # if mpi_rank == 0:
    #     file_desc, fname = tempfile.mkstemp(suffix=".dat", prefix="collect_bunch_", dir="/tmp")
    #     os.close(file_desc)
    #
    # TODO: this doesn't seem to work. "SystemError: PY_SSIZE_T_CLEAN macro must be defined for '#' formats"
    # fname = orbit_mpi.MPI_Bcast(fname, orbit_mpi.mpi_datatype.MPI_CHAR, 0, mpi_comm)

    # Using a fixed filename in the temp directory for now. Maybe that's sufficient.
    fname = f"/tmp/collect_bunch_tmpfile_{mpi_rank}.dat"

    local_shape = (local_rows, coords_shape[1])
    dtype = np.float64
    coords_memmap = np.memmap(fname, dtype=dtype, mode="w+", shape=local_shape)

    for i in range(local_rows):
        coords_memmap[i, BunchCoord.X] = bunch.x(i)
        coords_memmap[i, BunchCoord.XP] = bunch.xp(i)
        coords_memmap[i, BunchCoord.Y] = bunch.y(i)
        coords_memmap[i, BunchCoord.YP] = bunch.yp(i)
        coords_memmap[i, BunchCoord.Z] = bunch.z(i)
        coords_memmap[i, BunchCoord.DE] = bunch.dE(i)

    coords_memmap.flush()

    bunch_dict = {"coords": None, "sync_part": {}, "attributes": {}}

    if mpi_rank == 0:
        sync_part = bunch.getSyncParticle()

        bunch_dict["sync_part"] |= {
            "coords": np.array(sync_part.pVector()),
            "kin_energy": np.float64(sync_part.kinEnergy()),
            "momentum": np.float64(sync_part.momentum()),
            "beta": np.float64(sync_part.beta()),
            "gamma": np.float64(sync_part.gamma()),
            "time": np.float64(sync_part.time()),
        }

        for attr in bunch.bunchAttrDoubleNames():
            bunch_dict["attributes"][attr] = np.float64(bunch.bunchAttrDouble(attr))

        for attr in bunch.bunchAttrIntNames():
            bunch_dict["attributes"][attr] = np.int32(bunch.bunchAttrInt(attr))

    orbit_mpi.MPI_Barrier(mpi_comm)

    if mpi_rank == 0:
        coords_memmap = np.memmap(fname, dtype=dtype, mode="r+", shape=coords_shape)

        start_row = local_rows

        for r in range(1, orbit_mpi.MPI_Comm_size(mpi_comm)):
            src_fname = f"/tmp/collect_bunch_tmpfile_{r}.dat"

            if not os.path.exists(src_fname):
                raise FileNotFoundError(f"Expected temporary file '{src_fname}' not found. Something went wrong.")

            src_memmap = np.memmap(src_fname, dtype=dtype, mode="r")
            src_memmap = src_memmap.reshape((-1, coords_shape[1]))

            stop_row = start_row + src_memmap.shape[0]

            coords_memmap[start_row:stop_row, :] = src_memmap[:, :]
            coords_memmap.flush()

            del src_memmap
            os.remove(src_fname)
            start_row = stop_row

        bunch_dict["coords"] = (
            coords_memmap if return_memmap else np.array(coords_memmap)
        )

        return bunch_dict
