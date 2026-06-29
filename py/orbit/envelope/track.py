import numpy as np
import warnings

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle
from orbit.lattice import AccNode
from orbit.lattice import AccLattice
from orbit.teapot import BendTEAPOT
from orbit.py_linac.lattice import Bend as BendLINAC

from .matrix import get_matrix
from .envelope import Envelope


ENTRANCE = AccNode.ENTRANCE
BODY = AccNode.BODY
EXIT = AccNode.EXIT

BEFORE = AccNode.BEFORE
AFTER = AccNode.AFTER


class EnvelopeTracker:
    def __init__(self, lattice: AccLattice, sc: str | None = None) -> None:
        """Constructor.

        Args:
            lattice: The accelerator lattice.
            sc: Envelope space charge model {"2d", "3d", None}.
        """
        self.lattice = lattice
        self.sc = sc

        # For pre-computing elements
        self.elements = []
        self.one_turn_matrix = None
        # [TO DO] option to return one-turn matrix including linear space charge

        for node in self.lattice.getNodes():
            if type(node) in (BendTEAPOT, BendLINAC):
                if node.getParam("ea1") != 0.0 or node.getParam("ea2") != 0.0:
                    message = f"Found bend ea1 or ea2 != 0.0 ({node.getName()}.)"
                    message += " Nonzero edge angles are not yet supported in envelope tracking."
                    message += " Setting ea1 and ea2 to 0.0."
                    warnings.warn(message)

                    node.setParam("ea1", 0.0)
                    node.setParam("ea2", 0.0)

    def track(self, envelope: Envelope) -> None:
        sync_part = envelope.sync_part
        charge = envelope.charge()

        for node_index, node in enumerate(self.lattice.getNodes()):
            for child_node in node.getChildNodes(ENTRANCE):
                matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    envelope.transform(matrix)

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=BEFORE):
                    matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        envelope.transform(matrix)

                matrix_sc = None
                if self.sc:
                    length = node.getLength(part_index)
                    if length > 0:
                        if self.sc == "2d":
                            matrix_sc = envelope.sc_matrix_2d(length)
                        elif self.sc == "3d":
                            matrix_sc = envelope.sc_matrix_3d(length)
                        else:
                            raise ValueError

                matrix = get_matrix(node, sync_part=sync_part, charge=charge, index=part_index)
                if matrix is not None:
                    if matrix_sc is not None:
                        matrix = matrix @ matrix_sc
                    envelope.transform(matrix)

                for child_node in node.getChildNodes(BODY, part_index, place_in_part=AFTER):
                    matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        envelope.transform(matrix)

            for child_node in node.getChildNodes(EXIT):
                matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    envelope.transform(matrix)

    def track_history(self, envelope: Envelope) -> dict[str, list]:
        """Track and return envelope parameters vs. position in lattice."""
        history = {}
        history["position"] = []
        history["rms_x"] = []
        history["rms_y"] = []
        history["rms_z"] = []
        history["kin_energy"] = []

        sync_part = envelope.sync_part
        charge = envelope.charge()
        node_positions = self.lattice.getNodePositionsDict()

        history["position"].append(0.0)
        history["rms_x"].append(1000.0 * envelope.rms(0))
        history["rms_y"].append(1000.0 * envelope.rms(2))
        history["rms_z"].append(1000.0 * envelope.rms(4))
        history["kin_energy"].append(envelope.sync_part.kinEnergy())

        for node_index, node in enumerate(self.lattice.getNodes()):
            for child_node in node.getChildNodes(ENTRANCE):
                matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    envelope.transform(matrix)

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=BEFORE):
                    matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        envelope.transform(matrix)

                matrix_sc = None
                if self.sc:
                    length = node.getLength(part_index)
                    if length > 0:
                        if self.sc == "2d":
                            matrix_sc = envelope.sc_matrix_2d(length)
                        elif self.sc == "3d":
                            matrix_sc = envelope.sc_matrix_3d(length)
                        else:
                            raise ValueError

                matrix = get_matrix(node, sync_part=sync_part, charge=charge, index=part_index)
                if matrix is not None:
                    if matrix_sc is not None:
                        matrix = matrix @ matrix_sc
                    envelope.transform(matrix)

                position_start, position_stop = node_positions[node]
                position = position_start + node.getLength(part_index) * (part_index + 1)

                history["position"].append(position)
                history["rms_x"].append(1000.0 * envelope.rms(0))
                history["rms_y"].append(1000.0 * envelope.rms(2))
                history["rms_z"].append(1000.0 * envelope.rms(4))
                history["kin_energy"].append(envelope.sync_part.kinEnergy())

                for child_node in node.getChildNodes(BODY, part_index, place_in_part=AFTER):
                    matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        envelope.transform(matrix)

            for child_node in node.getChildNodes(EXIT):
                matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    envelope.transform(matrix)

        return history

    def precompute_matrices(self, envelope: Envelope) -> None:
        """Pre-compute transfer matrices for each node.

        For each node, return tuple (node, matrix). Mark space charge kicks as ("sc", length).
        """
        sync_part = envelope.sync_part
        charge = envelope.charge()

        self.elements = []
        for node_index, node in enumerate(self.lattice.getNodes()):
            for child_node in node.getChildNodes(ENTRANCE):
                matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    self.elements.append((child_node, matrix))

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(BODY, part_index, place_in_part=BEFORE):
                    matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        self.elements.append((child_node, matrix))

                if self.sc:
                    length = node.getLength(part_index)
                    if length > 0:
                        self.elements.append(("sc", length))

                matrix = get_matrix(node, sync_part=sync_part, charge=charge, index=part_index)
                if matrix is not None:
                    self.elements.append((node, matrix))

                for child_node in node.getChildNodes(BODY, part_index, place_in_part=AFTER):
                    matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                    if matrix is not None:
                        self.elements.append((node, matrix))

            for child_node in node.getChildNodes(EXIT):
                matrix = get_matrix(child_node, sync_part=sync_part, charge=charge)
                if matrix is not None:
                    self.elements.append((node, matrix))

    def track_ring(self, envelope: Envelope) -> None:
        """Track using pre-computed transfer matrices.

        The method assumes that all nodes are static and that there is no
        change in the synchronous particle energy. In this case the matrices
        can be computed once and reused on each turn. If there is no space charge,
        we track using the one-turn matrix.
        """

        # Pre-compute transfer matrices on the first turn.
        if not self.elements:
            self.precompute_matrices(envelope)
            self.one_turn_matrix = None

        # If there is no space charge, apply the one-turn transfer matrix.
        if not self.sc:
            if self.one_turn_matrix is None:
                self.one_turn_matrix = np.identity(7)
                for (node, matrix) in self.elements:
                    self.one_turn_matrix = np.matmul(self.one_turn_matrix, matrix)
            return envelope.transform(self.one_turn_matrix)

        # If there is space charge, apply the matrices one-by-one.
        for element in self.elements:
            if element[0] == "sc":
                length = element[1]
                if self.sc == "2d":
                    envelope.transform(envelope.sc_matrix_2d(length))
                elif self.sc == "3d":
                    envelope.transform(envelope.sc_matrix_3d(length))
                else:
                    raise ValueError
            else:
                node, matrix = element
                envelope.transform(matrix)
