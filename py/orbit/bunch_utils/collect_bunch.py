import os
import pathlib

from orbit.core.bunch import Bunch
from orbit.core import orbit_mpi
from orbit.bunch_utils.file_handler import FileHandler, NumPyHandler, BunchDict

import numpy as np


def collect_bunch(
    bunch: Bunch, output_dir: str | pathlib.Path = "/tmp", return_memmap: bool = True
) -> BunchDict | None:
    """Collects attributes from a PyOrbit Bunch across all MPI ranks and returns it as a dictionary.
    Parameters
    ----------
    bunch : Bunch
        The PyOrbit::Bunch object from which to collect attributes.
    output_dir : str | pathlib.Path, optional
        The director to use for temporary storage of the bunch coordinates on each MPI rank.
        If None, the bunch will be stored in "/tmp".
        Note: take care that the temporary files are created in a directory where all MPI ranks have write access.
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
                    and the 6 columns correspond to [x, xp, y, yp, z, dE] in units of [m, rad, m, rad, m, GeV], respectively.
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
    Raises
    ------
    FileNotFoundError
        If the temporary files created by non-root MPI ranks could not be found by the root rank during
        the collection process.
    """

    global_size = bunch.getSizeGlobal()

    if global_size == 0:
        return None

    mpi_comm = bunch.getMPIComm()
    mpi_rank = orbit_mpi.MPI_Comm_rank(mpi_comm)

    coords_shape = (bunch.getSizeGlobal(), 6)

    local_rows = bunch.getSize()

    if isinstance(output_dir, str):
        output_dir = pathlib.Path(output_dir)

    fname = output_dir / f"collect_bunch_tmpfile_{mpi_rank}.dat"

    local_shape = (local_rows, coords_shape[1])
    dtype = np.float64
    coords_memmap = np.memmap(fname, dtype=dtype, mode="w+", shape=local_shape)

    for i in range(local_rows):
        coords_memmap[i, :] = (
            bunch.x(i),
            bunch.xp(i),
            bunch.y(i),
            bunch.yp(i),
            bunch.z(i),
            bunch.dE(i),
        )

    coords_memmap.flush()

    bunch_dict: BunchDict = {"coords": None, "sync_part": {}, "attributes": {}}

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

    if mpi_rank != 0:
        return None

    coords_memmap = np.memmap(fname, dtype=dtype, mode="r+", shape=coords_shape)

    start_row = local_rows

    for r in range(1, orbit_mpi.MPI_Comm_size(mpi_comm)):
        src_fname = output_dir / f"collect_bunch_tmpfile_{r}.dat"

        if not os.path.exists(src_fname):
            raise FileNotFoundError(
                f"Expected temporary file '{src_fname}' not found. Something went wrong."
            )

        src_memmap = np.memmap(src_fname, dtype=dtype, mode="r")
        src_memmap = src_memmap.reshape((-1, coords_shape[1]))

        stop_row = start_row + src_memmap.shape[0]

        coords_memmap[start_row:stop_row, :] = src_memmap[:, :]
        coords_memmap.flush()

        del src_memmap
        os.remove(src_fname)
        start_row = stop_row

    bunch_dict["coords"] = coords_memmap if return_memmap else np.array(coords_memmap)

    return bunch_dict


def save_bunch(
    bunch: Bunch | BunchDict,
    output_dir: str | pathlib.Path = "bunch_data/",
    Handler: type[FileHandler] = NumPyHandler,
) -> None:
    """Saves the collected bunch attributes to a specified directory.
    Parameters
    ----------
    bunch_dict : Bunch | BunchDict
        The PyOrbit::Bunch object or the dictionary containing the collected bunch attributes.
    output_dir : str, optional
        The directory where the bunch data files will be saved. Default is "bunch_data/".
    Handler : FileHandler, optional
        The file handler class to use for writing the bunch data. Default is NumPyHandler.
    Returns
    -------
    None
    Raises
    ------
    ValueError
        If the provided `bunch` is neither a Bunch instance nor a BunchDict.
    """

    if isinstance(bunch, Bunch):
        mpi_comm = bunch.getMPIComm()
        bunch = collect_bunch(bunch)
    else:
        mpi_comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD

    mpi_rank = orbit_mpi.MPI_Comm_rank(mpi_comm)

    if mpi_rank != 0 or bunch is None:
        return

    if bunch["coords"].shape[0] == 0:
        print("No particles in the bunch to save.")
        return

    if isinstance(output_dir, str):
        output_dir = pathlib.Path(output_dir)

    handler = Handler(output_dir)
    handler.write(bunch)


def load_bunch(
    input_dir: str | pathlib.Path, Handler: type[FileHandler] = NumPyHandler
) -> tuple[Bunch, BunchDict]:
    """Loads the bunch attributes from a specified directory containing NumPy binary files.
    Parameters
    ----------
    input_dir : str | pathlib.Path
        The directory from which to load the bunch data files.
    Handler : FileHandler, optional
        The file handler class to use for reading the bunch data. Default is NumPyHandler.
        See `orbit.bunch_utils.file_handler` for available handlers.
    Returns
    -------
    BunchDict
        A dictionary containing the loaded bunch attributes.
    Raises
    ------
    FileNotFoundError
        If the required files are not found in the specified directory.
    TypeError
        If an attribute in the loaded bunch has an unsupported type.
    """
    mpi_comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
    mpi_rank = orbit_mpi.MPI_Comm_rank(mpi_comm)
    mpi_size = orbit_mpi.MPI_Comm_size(mpi_comm)

    handler = Handler(input_dir)

    bunch_dict = handler.read()

    coords = bunch_dict["coords"]

    global_size = coords.shape[0]

    local_size = global_size // mpi_size
    remainder = global_size % mpi_size
    if mpi_rank < remainder:
        local_size += 1
        start_row = mpi_rank * local_size
    else:
        start_row = mpi_rank * local_size + remainder
    stop_row = start_row + local_size

    local_coords = coords[start_row:stop_row, :]

    bunch = Bunch()

    for i in range(local_size):
        bunch.addParticle(*local_coords[i, :])

    for attr, value in bunch_dict["attributes"].items():
        if np.issubdtype(value, np.floating):
            bunch.bunchAttrDouble(attr, value)
        elif np.issubdtype(value, np.integer):
            bunch.bunchAttrInt(attr, value)
        else:
            raise TypeError(f"Unsupported attribute type for '{attr}': {type(value)}")

    sync_part_obj = bunch.getSyncParticle()
    sync_part_obj.rVector(tuple(bunch_dict["sync_part"]["coords"]))
    sync_part_obj.kinEnergy(bunch_dict["sync_part"]["kin_energy"])
    sync_part_obj.time(bunch_dict["sync_part"]["time"])

    return bunch, bunch_dict
