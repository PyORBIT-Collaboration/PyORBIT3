import pathlib
from typing import Any, Protocol, TypedDict

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
    coords: NDArray[np.float64]
    sync_part: SyncPartDict
    attributes: dict[str, np.float64 | np.int32]


class FileHandler(Protocol):
    """Protocol for file handlers to read/write bunch data."""

    def __init__(self, *args: Any, **kwargs: Any) -> None: ...

    def read(self) -> BunchDict: ...

    def write(self, bunch: BunchDict) -> None: ...


class NumPyHandler:
    """Handler implementing the FileHandler protocol for NumPy binary files.
    This handler will create two files in the directory passed to the constructor:
        - coords.npy: A memory-mapped NumPy array containing the bunch coordinates.
        - attributes.npz: A NumPy archive containing data related to the synchronous particle and other bunch attributes.
    """

    _coords_fname = "coords.npy"
    _attributes_fname = "attributes.npz"

    def __init__(self, dir_name: str | pathlib.Path):
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


class HDF5Handler:
    # TODO
    def __init__(self):
        raise NotImplementedError("HDF5Handler is not yet implemented.")

    def read(self) -> BunchDict:
        raise NotImplementedError("HDF5Handler is not yet implemented.")

    def write(self, bunch: BunchDict) -> None:
        raise NotImplementedError("HDF5Handler is not yet implemented.")
