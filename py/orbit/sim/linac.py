import math
import os
import sys
import time
from typing import Callable
from typing import Optional

import numpy as np

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.orbit_utils import BunchExtremaCalculator
from orbit.lattice import AccActionsContainer
from orbit.lattice import AccNode
from orbit.lattice import AccLattice
from orbit.utils.consts import speed_of_light


def get_z_to_phase_coeff(bunch: Bunch, frequency: float) -> float:
    wavelength = speed_of_light / frequency
    return -360.0 / (bunch.getSyncParticle().beta() * wavelength)


def reverse_bunch(bunch: Bunch) -> Bunch:
    size = bunch.getSize()
    for i in range(size):
        bunch.xp(i, -bunch.xp(i))
        bunch.yp(i, -bunch.yp(i))
        bunch.z(i, -bunch.z(i))
    return bunch
    

def track_bunch(
    bunch: Bunch,
    lattice: AccLattice,
    index_start: int = None,
    index_stop: int = None,
    copy: bool = False,
    **kwargs
) -> Bunch:
    """Track bunch forward or backward through the lattice."""
    if index_start is None:
        index_start = 0

    if index_stop is None:
        index_stop = len(lattice.getNodes()) - 1

    reverse = index_start > index_stop
    node_start = lattice.getNodes()[index_start]
    node_stop = lattice.getNodes()[index_stop]

    bunch_out = None
    if copy:
        bunch_out = Bunch()
        bunch.copyBunchTo(bunch_out)
    else:
        bunch_out = bunch

    if reverse:
        bunch_out = reverse_bunch(bunch_out)
        lattice.reverseOrder()

    lattice.trackBunch(
        bunch_out,
        index_start=lattice.getNodeIndex(node_start),
        index_stop=lattice.getNodeIndex(node_stop),
        **kwargs
    )

    if reverse:
        bunch_out = reverse_bunch(bunch_out)
        lattice.reverseOrder()

    return bunch_out



class BunchWriter:
    """Writes bunch to file.

    File name is "{output_dir}/bunch_{index}_{node_name}.dat".
    Example:
        - bunch_0001_QH05.dat
        - bunch_0002_QV06.dat
    """

    def __init__(self, output_dir: str = None, index: int = 0, verbose: int = 1) -> None:
        self.output_dir = output_dir
        self.index = index
        self.verbose = verbose
        self.position = 0.0

        self.mpi_comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        self.mpi_rank = orbit_mpi.MPI_Comm_rank(self.mpi_comm)

    def __call__(
        self, bunch: Bunch, node_name: str = None, position: float = None, filename: str = None
    ) -> None:
        if filename is None:
            filename = "bunch"
            if self.index is not None:
                filename = "{}_{:04.0f}".format(filename, self.index)
            if node_name is not None:
                node_name = node_name.replace(" ", "_")
                filename = "{}_{}".format(filename, node_name)
            filename = "{}.dat".format(filename)

        filename = os.path.join(self.output_dir, filename)

        if self.mpi_rank == 0 and self.verbose:
            print("Writing bunch to file {}".format(filename))

        bunch.dumpBunch(filename)

        if self.index is not None:
            self.index += 1

        if position is not None:
            self.position = position


class BunchMonitor:
    """Monitors bunch within linac."""

    def __init__(
        self,
        output_dir: str = None,
        stride: float = 0.1,
        stride_write: float = math.inf,
        position_offset: float = 0.0,
        rf_frequency: float = None,
        stop_node: Optional[str] = None,
        bunch_writer: BunchWriter = None,
        verbose: bool = True,
    ) -> None:
        """Constructor.

        Args:
            output_dir: Path to output directory.
            stride: Distance between scalar bunch measurements.
            stride_write: Distance between saving bunch to file.
            position_offset: Starting position in lattice [m].
            rf_frequency: For converting longitudinal position to phase.
            stop_node: Stop at this node if provided.
            bunch_writer: Writes bunch to file.
            verbose: Whether to print updates.
        """

        # Save MPI rank
        self.mpi_comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        self.mpi_rank = orbit_mpi.MPI_Comm_rank(self.mpi_comm)

        # Settings
        self.output_dir = output_dir
        self.stride = stride
        self.stride_write = stride_write
        self.rf_frequency = rf_frequency
        self.verbose = verbose
        self.stop_node = stop_node

        # State
        self.position = self.position_offset = position_offset
        self.index = 0
        self.start_time = None
        self.time_ellapsed = 0.0
        self.reached_stop_node = False

        # Helpers
        self.bunch_writer = bunch_writer

        # Store scalar history in `history` dictionary.
        if self.mpi_rank == 0:
            keys = [
                "position",
                "n_parts",
                "gamma",
                "beta",
                "energy",
                "x_rms",
                "y_rms",
                "z_rms",
                "z_rms_deg",
                "z_to_phase_coeff",
                "x_min",
                "x_max",
                "y_min",
                "y_max",
                "z_min",
                "z_max",
                "eps_x",
                "eps_y",
                "eps_z",
                "eps_xy",
                "eps_xyz",
            ]
            for i in range(6):
                keys.append("mean_{}".format(i))
            for i in range(6):
                for j in range(i + 1):
                    keys.append("cov_{}-{}".format(j, i))

            self.history = {}
            for key in keys:
                self.history[key] = []

            if self.rf_frequency is None:
                self.history.pop("z_rms_deg")
                self.history.pop("z_to_phase_coeff")

            if self.output_dir is not None:
                filename = os.path.join(self.output_dir, "history.dat")
                self.history_file = open(filename, "w")

                # Write header line
                header = "#"
                header = header + ",".join(keys)
                header = header[:-1] + "\n"
                self.history_file.write(header)

    def __call__(self, params_dict: dict, force_update: bool = False) -> None:
        """Measure the bunch.

        Args:
            params_dict: Dictionary with the following keys:
                "bunch": Reference to tracked Bunch object.
                "path_length": Total tracking distance.
                "node": Reference to current AccNode object.
            force_update: Forces measurement update.
        """
        # Update position; decide whether to proceed.
        position = params_dict["path_length"] + self.position_offset
        is_stop_node = (self.stop_node is not None) and (
            params_dict["node"].getName() == self.stop_node
        )

        if force_update:
            pass
        elif is_stop_node:
            if self.reached_stop_node:
                return
            self.reached_stop_node = True
        elif self.index > 0:
            if (position - self.position) < self.stride:
                return
        self.position = position

        # Update ellapsed time.
        if self.start_time is None:
            self.start_time = time.time()
        self.time_ellapsed = time.time() - self.start_time

        # Collect bunch and node from parameter dictionary.
        bunch = params_dict["bunch"]
        node = params_dict["node"]

        # Measure scalars.
        beta = bunch.getSyncParticle().beta()
        gamma = bunch.getSyncParticle().gamma()
        bunch_size_global = bunch.getSizeGlobal()
        if self.mpi_rank == 0:
            self.history["position"].append(position)
            self.history["n_parts"].append(bunch_size_global)
            self.history["gamma"].append(gamma)
            self.history["beta"].append(beta)
            self.history["energy"].append(bunch.getSyncParticle().kinEnergy())

        # Measure bunch centroid and 6 x 6 covariance matrix.
        twiss_analysis = BunchTwissAnalysis()
        twiss_analysis.computeBunchMoments(bunch, 2, 0, 0)

        centroid = np.zeros(6)
        for i in range(6):
            centroid[i] = twiss_analysis.getAverage(i)

        cov_matrix = np.zeros((6, 6))
        for i in range(6):
            for j in range(i + 1):
                cov_matrix[i, j] = cov_matrix[j, i] = twiss_analysis.getCorrelation(j, i)

        if self.mpi_rank == 0:
            for i in range(6):
                key = "mean_{}".format(i)
                value = centroid[i]
                self.history[key].append(value)

        if self.mpi_rank == 0:
            for i in range(6):
                for j in range(i + 1):
                    key = "cov_{}-{}".format(j, i)
                    value = cov_matrix[j, i]
                    self.history[key].append(value)

        # Record other parameters derived from covariance matrix.
        if self.mpi_rank == 0:
            x_rms = math.sqrt(cov_matrix[0, 0])
            y_rms = math.sqrt(cov_matrix[2, 2])
            z_rms = math.sqrt(cov_matrix[4, 4])
            self.history["x_rms"].append(x_rms)
            self.history["y_rms"].append(y_rms)
            self.history["z_rms"].append(z_rms)

            if self.rf_frequency is not None:
                z_to_phase_coeff = get_z_to_phase_coeff(bunch, self.rf_frequency)
                z_rms_deg = -z_to_phase_coeff * z_rms
                self.history["z_rms_deg"].append(z_rms_deg)
                self.history["z_to_phase_coeff"].append(z_to_phase_coeff)

            eps_x = np.sqrt(np.linalg.det(cov_matrix[0:2, 0:2]))
            eps_y = np.sqrt(np.linalg.det(cov_matrix[2:4, 2:4]))
            eps_z = np.sqrt(np.linalg.det(cov_matrix[4:6, 4:6]))
            eps_xy = np.sqrt(np.linalg.det(cov_matrix[0:4, 0:4]))
            eps_xyz = np.sqrt(np.linalg.det(cov_matrix))
            self.history["eps_x"].append(eps_x)
            self.history["eps_y"].append(eps_y)
            self.history["eps_z"].append(eps_z)
            self.history["eps_xy"].append(eps_xy)
            self.history["eps_xyz"].append(eps_xyz)

        # Measure min/max particle coordinates.
        extrema_calculator = BunchExtremaCalculator()
        (x_min, x_max, y_min, y_max, z_min, z_max) = extrema_calculator.extremaXYZ(bunch)
        if self.mpi_rank == 0:
            self.history["x_min"].append(x_min)
            self.history["x_max"].append(x_max)
            self.history["y_min"].append(y_min)
            self.history["y_max"].append(y_max)
            self.history["z_min"].append(z_min)
            self.history["z_max"].append(z_max)

        # Print update
        if self.verbose and (self.mpi_rank == 0):
            message = ""
            message += " index={:05.0f}".format(self.index)
            message += " t={:0.2f}".format(self.time_ellapsed)
            message += " s={:0.3f}".format(self.position)
            message += " xrms={:0.2f}".format(x_rms * 1000.0)
            message += " yrms={:0.2f}".format(y_rms * 1000.0)
            message += " zrms={:0.2f}".format(z_rms * 1000.0)
            message += " size={}".format(bunch_size_global)
            message += " node={}".format(node.getName())
            print(message)
            sys.stdout.flush()  # for MPI (bug?)

        # Increase index
        self.index += 1

        # Write phase space coordinates to file
        if self.bunch_writer is not None:
            if (position - self.bunch_writer.position) >= self.stride_write:
                self.bunch_writer(bunch, node_name=node.getName(), position=position)

        # Write new line to history file
        if (self.mpi_rank == 0) and (self.output_dir is not None):
            data = [self.history[key][-1] for key in self.history]
            line = ""
            for x in data:
                line = line + "{},".format(x)
            line = line[:-1] + "\n"
            self.history_file.write(line)
