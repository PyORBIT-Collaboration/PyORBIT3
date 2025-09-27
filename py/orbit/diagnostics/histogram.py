import os
import sys
import time
from typing import Any
from typing import Callable
from typing import Union

import numpy as np
import xarray as xr

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch
from orbit.core.spacecharge import Grid1D
from orbit.core.spacecharge import Grid2D
from orbit.core.spacecharge import Grid3D
from orbit.lattice import AccLattice
from orbit.lattice import AccNode


Grid = Union[Grid1D, Grid2D, Grid3D]


def get_grid_points(grid_coords: list[np.ndarray]) -> np.ndarray:
    if len(grid_coords) == 1:
        return grid_coords[0]
    return np.vstack([c.ravel() for c in np.meshgrid(*grid_coords, indexing="ij")]).T


def grid_edges_to_coords(grid_edges: np.ndarray) -> np.ndarray:
    return 0.5 * (grid_edges[:-1] + grid_edges[1:])


def make_grid(shape: tuple[int, ...], limits: list[tuple[float, float]]) -> Grid:

    ndim = len(shape)

    grid = None
    if ndim == 1:
        grid = Grid1D(shape[0] + 1, limits[0][0], limits[0][1])
    elif ndim == 2:
        grid = Grid2D(
            shape[0] + 1,
            shape[1] + 1,
            limits[0][0],
            limits[0][1],
            limits[1][0],
            limits[1][1],
        )
    elif ndim == 3:
        grid = Grid3D(
            shape[0] + 1,
            shape[1] + 1,
            shape[2] + 1,
        )
        grid.setGridX(limits[0][0], limits[0][1])
        grid.setGridY(limits[1][0], limits[1][1])
        grid.setGridZ(limits[2][0], limits[2][1])
    else:
        raise ValueError

    return grid


class BunchHistogram:
    """MPI-compatible bunch histogram."""
    def __init__(
        self,
        axis: tuple[int, ...],
        shape: tuple[int, ...],
        limits: list[tuple[float, float]],
        method: str = None,
        transform: Callable = None,
        normalize: bool = True,
        output_dir: str = None, 
        verbose: int = 2,
        **kwargs
    ) -> None:
        """Constructor.

        Args:
            axis: Axis on which to compute the histogram.
            shape: Number of bins along each axis.
            limits: Min/max coordinates along each axis.
            method: Smoothing method {"bilinear", "nine-point", None}.
            transform: Transforms bunch before histogram is calculated. 
                Call signature is `bunch_new = transform(bunch)`.
            normalize: Whehter to normalize values to PDF.
            output_dir: Output directory for saved files.
            verbose: Whether to print update messages.
        """
        self.mpi_comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        self.mpi_rank = orbit_mpi.MPI_Comm_rank(self.mpi_comm)
        self.output_dir = output_dir
        self.verbose = verbose

        self.axis = axis
        self.ndim = len(axis)
        self.method = method
        self.transform = transform
        self.normalize = normalize

        self.index = 0  # number of calls to `track` method
        self.node = None 

        if self.ndim > 2:
            raise NotImplementedError(
                "BunchHistogram does not yet support 3D grids. See "
                "https://github.com/PyORBIT-Collaboration/PyORBIT3/issues/46"
                " and "
                "https://github.com/PyORBIT-Collaboration/PyORBIT3/issues/47"
            )

        # Dimension names
        self.dims = ["x", "xp", "y", "yp", "z", "dE"]
        self.dims = [self.dims[i] for i in self.axis]

        # Create grid
        self.grid_shape = shape
        self.grid_limits = limits
        self.grid_edges = [
            np.linspace(self.grid_limits[i][0], self.grid_limits[i][1], self.grid_shape[i] + 1)
            for i in range(self.ndim)
        ]
        self.grid_coords = [grid_edges_to_coords(e) for e in self.grid_edges]
        self.grid_values = np.zeros(shape)
        self.grid_points = get_grid_points(self.grid_coords)
        self.grid = make_grid(self.grid_shape, self.grid_limits)

        # Store cell volume for normalization
        self.cell_volume = np.prod([e[1] - e[0] for e in self.grid_edges])

    def sync_mpi(self) -> None:
        self.grid.synchronizeMPI(self.mpi_comm)

    def bin_bunch(self, bunch: Bunch) -> None:
        macrosize = bunch.macroSize()
        if macrosize == 0:
            bunch.macroSize(1.0)

        if self.method == "bilinear":
            self.grid.binBunchBilinear(bunch, *self.axis)
        else:
            self.grid.binBunch(bunch, *self.axis)

        bunch.macroSize(macrosize)

    def compute_histogram(self, bunch: Bunch) -> np.ndarray:
        self.bin_bunch(bunch)
        self.sync_mpi()

        values = np.zeros(self.grid_points.shape[0])
        if self.method == "bilinear":
            for i, point in enumerate(self.grid_points):
                values[i] = self.grid.getValueBilinear(*point)
        elif self.method == "nine-point":
            for i, point in enumerate(self.grid_points):
                values[i] = self.grid.getValue(*point)
        else:
            for i, indices in enumerate(np.ndindex(*self.grid_shape)):
                values[i] = self.grid.getValueOnGrid(*indices)
        values = np.reshape(values, self.grid_shape)

        if self.normalize:
            values_sum = np.sum(values)
            if values_sum > 0.0:
                values /= values_sum
            values /= self.cell_volume
        return values

    def track(self, bunch: Bunch) -> None:
        bunch_copy = Bunch()
        bunch.copyBunchTo(bunch_copy)    
        if self.transform is not None:
            bunch_copy = self.transform(bunch_copy)

        self.grid.setZero()
        self.grid_values = self.compute_histogram(bunch_copy)

        if self.output_dir is not None:
            array = xr.DataArray(self.grid_values, coords=self.grid_coords, dims=self.dims)
            array.to_netcdf(path=self.get_filename())

        self.index += 1

    def get_filename(self) -> str:
        filename = "hist_" + "-".join([str(i) for i in self.axis])
        filename = "{}_{:04.0f}".format(filename, self.index)
        filename = "{}.nc".format(filename)
        filename = os.path.join(self.output_dir, filename)
        return filename


class BunchHistogram1D(BunchHistogram):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)


class BunchHistogram2D(BunchHistogram):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)


class BunchHistogram3D(BunchHistogram):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
