from __future__ import annotations

from dataclasses import dataclass
from typing import Any
from typing import TYPE_CHECKING

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle
from orbit.core.teapot_base import MatrixGenerator

from .envelope import Envelope
from .matrix_fitting import bunch_from_sync_particle
from .matrix_fitting import copy_sync_particle
from .matrix_fitting import fit_node_transfer_matrix

if TYPE_CHECKING:
    from orbit.lattice import AccActionsContainer
    from orbit.lattice import AccLattice
    from orbit.lattice import AccNode


ParamsDict = dict[str, Any]


@dataclass
class NodeMap:
    matrix: np.ndarray | None
    reference_position: tuple[float, float, float]
    reference_momentum: tuple[float, float, float]
    reference_x_axis: tuple[float, float, float]
    time_advance: float

    def update_sync_particle(self, sync_part: SyncParticle) -> None:
        current_time = sync_part.time()
        sync_part.rVector(self.reference_position)
        sync_part.pVector(self.reference_momentum)
        sync_part.nxVector(self.reference_x_axis)
        sync_part.time(current_time + self.time_advance)


class EnvelopeTracker:
    def __init__(
        self,
        lattice: AccLattice,
        index_start: int = 0,
        index_stop: int | None = None,
        fit: bool = True,
        fallback: bool = True,
    ) -> None:
        self.lattice = lattice

        self.index_start = index_start
        self.index_start = max(self.index_start, 0)

        self.index_stop = index_stop
        if self.index_stop is None or self.index_stop < 0:
            self.index_stop = len(self.lattice.getNodes())
        else:
            self.index_stop = min(self.index_stop + 1, len(self.lattice.getNodes()))

        self.fit = fit
        self.fallback = fallback
        
        self._node_maps = None
        self._reference_bunch = None
        self._matrix_generator = None

    def track(
        self,
        envelope: Envelope,
        sc: str | None = None,
        params_dict: ParamsDict | None = None,
        action_container: AccActionsContainer | None = None,
        return_matrix: bool = False
    ) -> np.ndarray | None:
        return self._run(
            envelope=envelope,
            sc=sc,
            params_dict=params_dict,
            action_container=action_container,
            return_matrix=return_matrix,
        )

    def precompute(self, envelope: Envelope, params_dict: ParamsDict | None = None) -> None:
        self._run(
            envelope=envelope.copy(),
            params_dict=None if params_dict is None else params_dict.copy(),
            precompute=True,
        )

    def _get_node_map(self, params_dict: ParamsDict) -> NodeMap:
        node = params_dict["node"]
        sync_part = params_dict["envelope"].sync_part
        if self._reference_bunch is None:
            self._reference_bunch = bunch_from_sync_particle(sync_part)
        reference_sync_part = self._reference_bunch.getSyncParticle()
        copy_sync_particle(sync_part, reference_sync_part)
        start_time = sync_part.time()
        part_index = node.getActivePartIndex()
        use_fit = self.fit

        if not use_fit:
            try:
                matrix = node.getMatrix(reference_sync_part, part_index=part_index)
            except NotImplementedError:
                if not self.fallback:
                    raise
                copy_sync_particle(sync_part, reference_sync_part)
                use_fit = True

        if use_fit:
            if self._matrix_generator is None:
                self._matrix_generator = MatrixGenerator()
            matrix = fit_node_transfer_matrix(
                node=node,
                bunch=self._reference_bunch,
                part_index=part_index,
                parent_node=params_dict["parentNode"],
                matrix_generator=self._matrix_generator,
                params_dict=params_dict,
            )

        return NodeMap(
            matrix=matrix,
            reference_position=tuple(reference_sync_part.rVector()),
            reference_momentum=tuple(reference_sync_part.pVector()),
            reference_x_axis=tuple(reference_sync_part.nxVector()),
            time_advance=reference_sync_part.time() - start_time,
        )

    def _run(
        self,
        envelope: Envelope,
        sc: str | None = None,
        params_dict: ParamsDict | None = None,
        action_container: AccActionsContainer | None = None,
        return_matrix: bool = False,
        precompute: bool = False,
    ) -> np.ndarray | None:

        from orbit.lattice import AccActionsContainer
        from orbit.lattice import AccNode

        if sc not in (None, False, "2d", "3d"):
            raise ValueError(f"Invalid envelope space charge option `{sc}`")

        if params_dict is None:
            params_dict = {}

        if action_container is None:
            action_container = AccActionsContainer()

        params_dict.update(
            envelope=envelope,
            lattice=self.lattice,
            action_container=action_container
        )
        params_dict.setdefault("path_length", 0.0)
        if hasattr(self.lattice, "getUseRealCharge"):
            params_dict.setdefault("useCharge", self.lattice.getUseRealCharge())

        total_matrix = None
        if return_matrix:
            total_matrix = np.eye(7)

        self._reference_bunch = None
        self._matrix_generator = None

        node_maps = {} if precompute else self._node_maps

        def apply_matrix(matrix: np.ndarray | None) -> None:
            nonlocal total_matrix
            if matrix is not None:
                envelope.transform(matrix)
                if total_matrix is not None:
                    total_matrix = matrix @ total_matrix

        def push_envelope(params_dict: ParamsDict) -> None:
            node = params_dict["node"]
            length = node.getLength(node.getActivePartIndex())
            if sc and length > 0:
                apply_matrix(envelope.sc_matrix(length, model=sc))
            part_index = node.getActivePartIndex()
            if node_maps is not None and not precompute:
                node_map = node_maps[node, part_index]
            else:
                node_map = self._get_node_map(params_dict)
                if precompute:
                    node_maps[node, part_index] = node_map
            node_map.update_sync_particle(envelope.sync_part)
            apply_matrix(node_map.matrix)

        action_container.addAction(push_envelope, AccNode.BODY)
        if self.index_start < self.index_stop:
            self.lattice.trackActions(
                action_container,
                params_dict,
                self.index_start,
                self.index_stop - 1,
            )
        action_container.removeAction(push_envelope, AccNode.BODY)
        if precompute:
            self._node_maps = node_maps
        return total_matrix
