from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle
from orbit.core.orbit_utils import Matrix
from orbit.core.teapot_base import MatrixGenerator

if TYPE_CHECKING:
    from orbit.lattice import AccLattice
    from orbit.lattice import AccNode


def orbit_matrix_to_numpy(matrix: Matrix) -> np.ndarray:
    matrix_out = np.zeros(matrix.size())
    for i in range(matrix_out.shape[0]):
        for j in range(matrix_out.shape[1]):
            matrix_out[i, j] = matrix.get(i, j)
    return matrix_out


def copy_sync_particle(source: SyncParticle, target: SyncParticle) -> None:
    target.rVector(source.rVector())
    target.pVector(source.pVector())
    target.nxVector(source.nxVector())
    target.time(source.time())


def bunch_from_sync_particle(sync_part: SyncParticle) -> Bunch:
    bunch = Bunch()
    bunch.mass(sync_part.mass())
    bunch.charge(sync_part.charge())
    copy_sync_particle(sync_part, bunch.getSyncParticle())
    return bunch


def fit_node_transfer_matrix(
    node: AccNode,
    bunch: Bunch,
    part_index: int = -1,
    parent_node: AccNode | AccLattice | None = None,
    matrix_generator: MatrixGenerator | None = None,
    lost_bunch: Bunch | None = None,
    params_dict: dict | None = None,
) -> np.ndarray:

    if matrix_generator is None:
        matrix_generator = MatrixGenerator()
    if lost_bunch is None:
        lost_bunch = Bunch()
        bunch.copyEmptyBunchTo(lost_bunch)

    matrix = Matrix(7, 7)
    matrix.unit()
    matrix_generator.initBunch(bunch)

    params = {} if params_dict is None else params_dict.copy()
    params["bunch"] = bunch
    params["lostbunch"] = lost_bunch
    params["node"] = node
    if parent_node is not None:
        params["parentNode"] = parent_node

    active_part_index = node.getActivePartIndex()
    if part_index >= 0:
        node.setActivePartIndex(part_index)
    try:
        node.track(params)
    finally:
        node.setActivePartIndex(active_part_index)

    matrix_generator.calculateMatrix(bunch, matrix)
    return orbit_matrix_to_numpy(matrix)


def fit_transfer_matrix(
    lattice: AccLattice,
    bunch: Bunch,
    index_start: int = 0,
    index_stop: int | None = None,
) -> np.ndarray:

    if index_stop is None:
        index_stop = -1

    matrix = Matrix(7, 7)
    matrix.unit()

    bunch_out = Bunch()
    lost_bunch = Bunch()
    bunch.copyEmptyBunchTo(bunch_out)
    bunch.copyEmptyBunchTo(lost_bunch)

    matrix_generator = MatrixGenerator()
    matrix_generator.initBunch(bunch_out)
    lattice.trackBunch(
        bunch_out,
        paramsDict={"lostbunch": lost_bunch},
        index_start=index_start,
        index_stop=index_stop,
    )
    matrix_generator.calculateMatrix(bunch_out, matrix)
    return orbit_matrix_to_numpy(matrix)
