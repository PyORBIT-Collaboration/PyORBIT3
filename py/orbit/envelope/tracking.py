from __future__ import annotations

from dataclasses import dataclass
from typing import Any
from typing import TYPE_CHECKING

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.teapot_base import MatrixGenerator

from .envelope import Envelope
from .node_matrices import get_node_matrix
from .matrix_fitting import fit_node_matrix

if TYPE_CHECKING:
    from orbit.lattice import AccActionsContainer
    from orbit.lattice import AccLattice
    from orbit.lattice import AccNode


ParamsDict = dict[str, Any]


@dataclass
class NodeMap:
    matrix: np.ndarray | None
    kin_energy: float
    time: float

    def update_envelope(self, envelope: Envelope) -> None:
        envelope.kin_energy = self.kin_energy
        envelope.time = self.time


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
        self.index_start = max(index_start, 0)
        self.index_stop = index_stop
        if self.index_stop is None or self.index_stop < 0:
            self.index_stop = len(self.lattice.getNodes())
        else:
            self.index_stop = min(self.index_stop + 1, len(self.lattice.getNodes()))

        self.fit = fit
        self.fallback = fallback
        
        self._bunch = None
        self._matrix_generator = None
        self._node_maps = None

    def track(
        self,
        envelope: Envelope,
        sc: str | None = None,
        params_dict: ParamsDict | None = None,
        action_container: AccActionsContainer | None = None,
        return_matrix: bool = False,
        store: bool = False,
    ) -> np.ndarray | None:

        from orbit.lattice import AccActionsContainer
        from orbit.lattice import AccNode

        if sc not in (None, False, "2d", "3d"):
            raise ValueError(f"Invalid envelope space charge option `{sc}`")

        if params_dict is None:
            params_dict = {}

        if action_container is None:
            action_container = AccActionsContainer()

        params_dict["lattice"] = self.lattice
        params_dict["envelope"] = envelope
        params_dict.setdefault("path_length", 0.0)
        if hasattr(self.lattice, "getUseRealCharge"):
            params_dict.setdefault("useCharge", self.lattice.getUseRealCharge())

        total_matrix = np.eye(7) if return_matrix else None
        node_maps = {} if store else self._node_maps
        identity_matrix = np.eye(7) if store else None

        def apply_matrix(matrix: np.ndarray | None) -> None:
            nonlocal total_matrix
            if matrix is not None:
                envelope.apply_matrix(matrix)
                if total_matrix is not None:
                    total_matrix = matrix @ total_matrix

        def push_envelope(params_dict: ParamsDict) -> None:
            node = params_dict["node"]
            part_index = node.getActivePartIndex()
            if sc:
                length = node.getLength(part_index)
                if length > 0:
                    apply_matrix(envelope.sc_matrix(length, model=sc))
            if node_maps and not store:
                node_map = node_maps[node, part_index]
                node_map.update_envelope(envelope)
            else:
                node_map = self._get_node_map(params_dict)
                if store:
                    if node_map.matrix is not None and np.array_equal(node_map.matrix, identity_matrix):
                        node_map.matrix = None
                    node_maps[node, part_index] = node_map
            apply_matrix(node_map.matrix)

        action_container.addAction(push_envelope, AccNode.BODY)
        try:
            if self.index_start < self.index_stop:
                self.lattice.trackActions(
                    action_container,
                    params_dict,
                    self.index_start,
                    self.index_stop - 1,
                )
        finally:
            action_container.removeAction(push_envelope, AccNode.BODY)
        if store:
            self._node_maps = node_maps
        return total_matrix

    def precompute(self, envelope: Envelope, params_dict: ParamsDict | None = None) -> None:
        self.track(
            envelope=envelope.copy(),
            params_dict=None if params_dict is None else params_dict.copy(),
            store=True,
        )

    def _get_node_map(self, params_dict: ParamsDict) -> NodeMap:
        node = params_dict["node"]
        envelope = params_dict["envelope"]
        initial_kinetic_energy = envelope.kin_energy
        initial_time = envelope.time

        part_index = node.getActivePartIndex()
        use_fit = self.fit

        if not use_fit:
            try:
                matrix = get_node_matrix(node, envelope, part_index=part_index)
            except NotImplementedError:
                envelope.kin_energy = initial_kinetic_energy
                envelope.time = initial_time
                if not self.fallback:
                    raise
                use_fit = True

        if use_fit:
            if self._bunch is None:
                self._bunch = Bunch()
            self._bunch.mass(envelope.mass)
            self._bunch.charge(envelope.charge)
            sync_particle = self._bunch.getSyncParticle()
            sync_particle.kinEnergy(envelope.kin_energy)
            sync_particle.time(envelope.time)
            if self._matrix_generator is None:
                self._matrix_generator = MatrixGenerator()
            matrix = fit_node_matrix(
                node=node,
                bunch=self._bunch,
                part_index=part_index,
                matrix_generator=self._matrix_generator,
                params_dict=params_dict,
            )
            envelope.kin_energy = sync_particle.kinEnergy()
            envelope.time = sync_particle.time()

        return NodeMap(
            matrix=matrix,
            kin_energy=envelope.kin_energy,
            time=envelope.time,
        )
