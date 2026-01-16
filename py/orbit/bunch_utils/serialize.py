import os
import pathlib
from typing import Any, Protocol, TypedDict

import numpy as np
from numpy.typing import NDArray

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch


class SyncPartDict(TypedDict):
    """A dictionary containing the attributes of the synchronous particle.

    Attributes
    ----------
    coords : NDArray[np.float64]
        The coordinates (x, px, y, py, z, dE) of the synchronous particle.
    kin_energy : float
        The kinetic energy of the synchronous particle.
    momentum : float
        The momentum of the synchronous particle.
    beta : float
        The beta of the synchronous particle.
    gamma : float
        The gamma of the synchronous particle.
    time : float
        The time of the synchronous particle.
    """
    coords: NDArray[np.float64]
    kin_energy: np.float64
    momentum: np.float64
    beta: np.float64
    gamma: np.float64
    time: np.float64


class BunchDict(TypedDict):
    """A dictionary containing the attributes of a PyOrbit::Bunch object.

    Attributes
    ----------
    coords : NDArray[np.float64]
        The coordinates (x, px, y, py, z, dE) of the particles in the bunch.
    sync_part : SyncPartDict
        The attributes of the synchronous particle.
    attributes : dict[str, np.float64 | np.int32]
        Other attributes of the bunch.
    """
    coords: NDArray[np.float64]
    sync_part: SyncPartDict
    attributes: dict[str, np.float64 | np.int32]


class FileHandler(Protocol):
    """Protocol for file handlers to read/write bunch data.

    Methods
    _______
    read() -> BunchDict:
        Reads the bunch data from the specified directory and returns it as a dictionary.
    write(bunch: BunchDict) -> None:
        Writes the bunch data to the specified directory.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None: ...

    def read(self) -> BunchDict: ...

    def write(self, bunch: BunchDict) -> None: ...


class NumPyHandler:
    """Handler implementing the FileHandler protocol for NumPy binary files.

    Attributes
    ----------
    _coords_fname : str
        The name of the file containing the bunch coordinates (default: "coords.npy").
    _attributes_fname : str
        The name of the file containing data related to the synchronous particle and other bunch attributes (default: "attributes.npz").
    """

    _coords_fname = "coords.npy"
    _attributes_fname = "attributes.npz"

    def __init__(self, dir_name: str | os.PathLike):
        """Constructor for the NumPyHandler class.

        Parameters
        ----------
        dir_name
            The directory in which to store the bunch data.

        """

        if isinstance(dir_name, str):
            dir_name = pathlib.Path(dir_name)
        self._dir_name = dir_name
        self._coords_path = dir_name / self._coords_fname
        self._attributes_path = dir_name / self._attributes_fname

    def read(self) -> BunchDict:
        if not self._coords_path.exists() or not self._attributes_path.exists():
            raise FileNotFoundError(
                f"Required files not found in directory: {self._dir_name}"
            )

        coords = np.load(self._coords_path, mmap_mode="r")

        attr_data = np.load(self._attributes_path, allow_pickle=True)

        sync_part = attr_data["sync_part"].item()
        attributes = attr_data["attributes"].item()

        return BunchDict(coords=coords, sync_part=sync_part, attributes=attributes)

    def write(self, bunch: BunchDict) -> None:
        self._dir_name.mkdir(parents=True, exist_ok=True)
        np.save(self._coords_path, bunch["coords"])
        np.savez(
            self._attributes_path,
            sync_part=bunch["sync_part"],
            attributes=bunch["attributes"],
        )


def collect_bunch(
    bunch: Bunch, output_dir: str | os.PathLike = "/tmp", return_memmap: bool = True
) -> BunchDict | None:
    """Collects attributes from a PyOrbit Bunch across all MPI ranks and returns it as a dictionary.

    Parameters
    ----------
    bunch : Bunch
        The PyOrbit::Bunch object from which to collect attributes.
    output_dir : str | os.PathLike, optional
        The director to use for temporary storage of the bunch coordinates on each MPI rank (default: "/tmp").
    return_memmap : bool, optional
        Return the bunch coordinates as a memory-mapped NumPy array, otherwise the
        entire array is copied into RAM and returned as normal NDArray (default: True).

    Note
    ----
    Take care that the temporary files are created in a directory which all MPI ranks have write access.

    Returns
    -------
    BunchDict | None
        A dictionary containing the collected bunch attributes. Returns None if not on the root MPI rank or if the global bunch size is 0.

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
    output_dir: str | os.PathLike = "bunch_data/",
    Handler: type[FileHandler] = NumPyHandler,
) -> None:
    """Saves the collected bunch attributes to a specified directory.

    Parameters
    ----------
    bunch_dict : Bunch | BunchDict
        The PyOrbit::Bunch object or the dictionary containing the collected bunch attributes.
    output_dir : str, optional
        The directory where the bunch data files will be saved (default: "bunch_data/").
    Handler : FileHandler, optional
        The file handler class to use for writing the bunch data (default: NumPyHandler).

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
    input_dir: str | os.PathLike, Handler: type[FileHandler] = NumPyHandler
) -> tuple[Bunch, BunchDict]:
    """Loads the bunch attributes from a specified directory containing NumPy binary files.

    Parameters
    ----------
    input_dir : str | os.PathLike
        The directory from which to load the bunch data files.
    Handler : FileHandler, optional
        The file handler class to use for reading the bunch data (default: NumPyHandler).
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
