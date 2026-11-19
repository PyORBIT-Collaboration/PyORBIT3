"""Compare HEBT1 phase advances from beta integration and transport matrices."""

import math
from pathlib import Path
import random

import pytest

from orbit.bunch_generators import TwissContainer, WaterBagDist3D
from orbit.core.bunch import Bunch, BunchTwissAnalysis
from orbit.core.orbit_mpi import (
    MPI_Bcast, MPI_Comm_rank, MPI_Comm_size, mpi_comm, mpi_datatype,
)
from orbit.lattice import AccActionsContainer
from orbit.py_linac.lattice import LinacAccLattice, LinacTrMatricesController
from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory
from orbit.utils import phaseNearTargetPhaseDeg


@pytest.fixture
def hebt1_lattice() -> LinacAccLattice:
    xml = Path(__file__).resolve().parent / "inputs" / "sns_hebt1.xml"
    factory = SNS_LinacLatticeFactory()
    factory.setMaxDriftLength(0.01)
    return factory.getLinacAccLattice(["HEBT1"], str(xml))


@pytest.fixture
def hebt1_bunch() -> Bunch:
    """Build a zero-current H- bunch.

    Note:
        A fixed seed makes particle sampling reproducible. The try-finally
        block restores the previous global RNG state, even if setup fails,
        so this fixture does not affect subsequent tests in the same process.
    """
    initial = Bunch()
    initial.mass(0.939294)
    initial.charge(-1.0)
    initial.getSyncParticle().kinEnergy(1.000222)
    initial.getSyncParticle().time(0.0)
    initial.macroSize(0.0)
    distribution = WaterBagDist3D(
        TwissContainer(1.4040, 5.6534, 0.1908e-6),
        TwissContainer(-2.3665, 13.0151, 0.1544e-6),
        TwissContainer(-16.436, 169.5695, 0.4265e-6),
    )
    comm = mpi_comm.MPI_COMM_WORLD
    rank, size = MPI_Comm_rank(comm), MPI_Comm_size(comm)
    random_state = random.getstate()
    try:
        random.seed(100)
        for index in range(1000):
            coordinates = distribution.getCoordinates()
            coordinates = MPI_Bcast(coordinates, mpi_datatype.MPI_DOUBLE, 0, comm)
            if index % size == rank:
                initial.addParticle(*coordinates)
    finally:
        random.setstate(random_state)
    return initial


def test_hebt1_phase_advance_from_transport_matrices(
    hebt1_lattice: LinacAccLattice, hebt1_bunch: Bunch
) -> None:
    """Check that beta-integrated and matrix-derived HEBT1 phase advances agree.

    Compare the entrance-to-exit phase advances in X, Y, and Z for a
    zero-current beam, requiring a Euclidean difference of at most 1 degree.
    """
    controller = LinacTrMatricesController()
    nodes = hebt1_lattice.getNodes()
    # The zero-length start/end markers span the complete HEBT1 sequence.
    matrices = controller.addTrMatrixGenNodes(hebt1_lattice, [nodes[0], nodes[-1]])
    for matrix in matrices:
        matrix.setTwissWeightUse(True, True, True)

    analysis = BunchTwissAnalysis()
    integrated_phase_advances = [0.0, 0.0, 0.0]
    params = {"old_pos": 0.0}
    actions = AccActionsContainer("HEBT1 phase integration")

    def _integrate(params):
        position = params["path_length"]
        delta = position - params["old_pos"]
        params["old_pos"] = position
        if delta == 0.0:
            return
        tracked = params["bunch"]
        sync = tracked.getSyncParticle()
        analysis.analyzeBunch(tracked)
        longitudinal_factor = 1.0 / (sync.beta() ** 2 * sync.gamma() ** 3 * tracked.mass())
        for axis in range(3):
            _, beta, _, _ = analysis.getTwiss(axis)
            factor = longitudinal_factor if axis == 2 else 1.0
            integrated_phase_advances[axis] = phaseNearTargetPhaseDeg(
                integrated_phase_advances[axis] + math.degrees(factor * delta / beta), 0.0
            )

    for location in (AccActionsContainer.ENTRANCE, AccActionsContainer.BODY, AccActionsContainer.EXIT):
        actions.addAction(_integrate, location)
    hebt1_lattice.trackBunch(hebt1_bunch, paramsDict=params, actionContainer=actions)

    matrix_phase_advances = matrices[-1].getPhaseAdvances()
    difference = math.dist(matrix_phase_advances, integrated_phase_advances)
    assert difference <= 1.0, (
        f"matrix phases {matrix_phase_advances}, "
        f"integrated phases {integrated_phase_advances}, "
        f"Euclidean difference {difference} degrees"
    )
