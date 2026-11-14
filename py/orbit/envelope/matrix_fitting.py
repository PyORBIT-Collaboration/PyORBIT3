from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.orbit_utils import Matrix
from orbit.core.teapot_base import MatrixGenerator
if TYPE_CHECKING:
    from orbit.lattice import AccLattice
    from orbit.lattice import AccNode


def orbit_matrix_to_numpy(matrix: Matrix) -> np.ndarray:
    matrix_out = np.empty(matrix.size(), dtype=np.float64)
    matrix.copyToBuffer(matrix_out)
    return matrix_out


def fit_node_matrix(
    node: AccNode,
    bunch: Bunch,
    part_index: int = -1,
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

    params_dict = {} if params_dict is None else params_dict.copy()
    params_dict["bunch"] = bunch
    params_dict["lostbunch"] = lost_bunch
    params_dict["node"] = node

    active_part_index = node.getActivePartIndex()
    if part_index >= 0:
        node.setActivePartIndex(part_index)
    try:
        node.track(params_dict)
    finally:
        node.setActivePartIndex(active_part_index)

    matrix_generator.calculateMatrix(bunch, matrix)
    return orbit_matrix_to_numpy(matrix)


def fit_lattice_matrix(
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
