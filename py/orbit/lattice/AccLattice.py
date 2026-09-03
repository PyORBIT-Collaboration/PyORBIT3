from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np

from ..utils import orbitFinalize
from ..utils import NamedObject
from ..utils import TypedObject

from .AccActionsContainer import AccActionsContainer
from .AccNode import AccNode

if TYPE_CHECKING:
    from orbit.envelope.envelope import Envelope


class AccLattice(NamedObject, TypedObject):
    """
    Class. The accelerator lattice class contains child nodes.
    """

    ENTRANCE = AccActionsContainer.ENTRANCE
    BODY = AccActionsContainer.BODY
    EXIT = AccActionsContainer.EXIT

    BEFORE = AccActionsContainer.BEFORE
    AFTER = AccActionsContainer.AFTER

    def __init__(self, name="no name"):
        """
        Constructor. Creates an empty accelerator lattice.
        """
        NamedObject.__init__(self, name)
        TypedObject.__init__(self, "lattice")
        self.__length = 0.0
        self.__isInitialized = False
        self.__children = []
        self.__childPositions = {}
        self.__envelopeElements = []
        self.__envelopeOneTurnMatrix = None
        self.__envelopeSpaceCharge = None

    def initialize(self):
        """
        Method. Initializes the lattice and child node structures.
        """
        res_dict = {}
        for node in self.__children:
            if node in res_dict:
                msg = "The AccLattice class instance should not have duplicate nodes!"
                msg = msg + os.linesep
                msg = msg + "Method initialize():"
                msg = msg + os.linesep
                msg = msg + "Name of node=" + node.getName()
                msg = msg + os.linesep
                msg = msg + "Type of node=" + node.getType()
                msg = msg + os.linesep
                orbitFinalize(msg)
            else:
                res_dict[node] = None
            node.initialize()
        del res_dict

        paramsDict = {}
        actions = AccActionsContainer()
        d = [0.0]
        posn = {}

        def accNodeExitAction(paramsDict):
            """
            Nonbound function. Sets lattice length and node
            positions. This is a closure (well, maybe not
            exactly). It uses external objects.
            """
            node = paramsDict["node"]
            parentNode = paramsDict["parentNode"]
            if isinstance(parentNode, AccLattice):
                posBefore = d[0]
                d[0] += node.getLength()
                posAfter = d[0]
                posn[node] = (posBefore, posAfter)

        actions.addAction(accNodeExitAction, AccNode.EXIT)
        self.trackActions(actions, paramsDict)
        self.__length = d[0]
        self.__childPositions = posn
        self.__isInitialized = True

    def isInitialized(self):
        """
        Method. Returns the initialization status (True or False).
        """
        return self.__isInitialized

    def addNode(self, node, index=-1):
        """
        Method. Adds a child node into the lattice. If the user
        specifies the index >= 0 the element will be inserted in
        the specified position into the children array
        """
        if isinstance(node, AccNode) == True:
            if index < 0:
                self.__children.append(node)
            else:
                self.__children.insert(index, node)
            self.__isInitialized = False

    def getNodes(self):
        """
        Method. Returns a list of all children
        of the first level in the lattice.
        """
        return self.__children

    def setNodes(self, childrenNodes):
        """
        Method. Set up a new list of all children
        of the first level in the lattice.
        """
        self.__children = childrenNodes

    def getNodeForName(self, name):
        """
        Method. Returns the node with certain name.
        """
        nodes = []
        for node in self.__children:
            if node.getName() == name:
                nodes.append(node)
        if len(nodes) == 1:
            return nodes[0]
        else:
            if len(nodes) == 0:
                return None
            else:
                msg = "The AccLattice class. Method getNodeForName found many nodes instead of one!"
                msg = msg + os.linesep
                msg = msg + "looking for name=" + name
                msg = msg + os.linesep
                msg = msg + "found nodes:"
                for node in nodes:
                    msg = msg + " " + node.getName()
                msg = msg + os.linesep
                msg = msg + "Please use getNodesForName method instead."
                msg = msg + os.linesep
                orbitFinalize(msg)

    def getNodesForName(self, name):
        """
        Method. Returns nodes with a certain name.
        """
        nodes = []
        for node in self.__children:
            if node.getName().find(name) == 0:
                nodes.append(node)
        return nodes

    def getNodesOfClass(self, class_of_node):
        """
        Method. Returns nodes off a certain class.
        """
        nodes = []
        for node in self.__children:
            if isinstance(node, class_of_node):
                nodes.append(node)
        return nodes

    def getNodesForSubstring(self, sub, no_sub=None):
        """
        Method. Returns nodes with names each of them has the certain substring.
        It is also possible to specify the unwanted substring as no_sub parameter.
        """
        nodes = []
        for node in self.__children:
            if no_sub == None:
                if node.getName().find(sub) >= 0:
                    nodes.append(node)
            else:
                if node.getName().find(sub) >= 0 and node.getName().find(no_sub) < 0:
                    nodes.append(node)
        return nodes

    def getNodeIndex(self, node):
        """
        Method. Returns the index of the node in the upper level of the lattice children-nodes.
        """
        return self.__children.index(node)

    def getNodePositionsDict(self):
        """
        Method. Returns a dictionary of
        {node:(start position, stop position)}
        tuples for all children of the first level in the lattice.
        """
        return self.__childPositions

    def getLength(self):
        """
        Method. Returns the physical length of the lattice.
        """
        return self.__length

    def reverseOrder(self):
        """
        This method is used for a lattice reversal and a bunch backtracking.
        This method will reverse the order of the children nodes. It will
        apply the reverse recursively to the all children nodes.
        """
        self.__children.reverse()
        for node in self.__children:
            node.reverseOrder()
        self.initialize()

    def structureToText(self):
        """
        Returns the text with the lattice structure.
        """
        txt = "==== START ==== Lattice =" + self.getName() + "  L=" + str(self.getLength())
        txt += os.linesep
        for node in self.__children:
            txt += node.structureToText("")
        txt += "==== STOP  ==== Lattice =" + self.getName() + "  L=" + str(self.getLength())
        txt += os.linesep
        return txt

    def _getSubLattice(self, accLatticeNew, index_start=-1, index_stop=-1):
        """
        It returns the sub-accelerator lattice with children with
        indexes between index_start and index_stop, inclusive. The
        subclasses of AccLattice should NOT override this method.
        """
        if index_start < 0:
            index_start = 0
        if index_stop < 0:
            index_stop = len(self.__children) - 1
        # clear the node array in the new sublattice
        accLatticeNew.setNodes([])
        for node in self.__children[index_start : index_stop + 1]:
            accLatticeNew.addNode(node)
        accLatticeNew.initialize()
        return accLatticeNew

    def getSubLattice(
        self,
        index_start=-1,
        index_stop=-1,
    ):
        """
        It returns the sub-accelerator lattice with children with
        indexes between index_start and index_stop inclusive. The
        subclasses of AccLattice should override this method to replace
        AccLattice() constructor by the sub-class type constructor
        """
        return self._getSubLattice(AccLattice(), index_start, index_stop)

    def trackActions(self, actionsContainer, paramsDict={}, index_start=-1, index_stop=-1):
        """
        Method. Tracks the actions through all nodes in the lattice. The indexes are inclusive.
        """
        paramsDict["lattice"] = self
        paramsDict["actions"] = actionsContainer
        if not ("path_length" in paramsDict):
            paramsDict["path_length"] = 0.0
        if index_start < 0:
            index_start = 0
        if index_stop < 0:
            index_stop = len(self.__children) - 1
        for node in self.__children[index_start : index_stop + 1]:
            paramsDict["node"] = node
            paramsDict["parentNode"] = self
            node.trackActions(actionsContainer, paramsDict)

    def _getNodesInRange(self, index_start: int = 0, index_stop: int = None):
        if index_stop is None:
            index_stop = len(self.__children) - 1
        return self.__children[index_start : index_stop + 1]

    def _prepareEnvelopeTracking(self) -> None:
        for node in self.__children:
            node_type = type(node)
            is_teapot_bend = node_type.__name__ == "BendTEAPOT" and node_type.__module__ == "orbit.teapot.teapot"
            is_linac_bend = node_type.__name__ == "Bend" and node_type.__module__ == "orbit.py_linac.lattice.LinacAccNodes"
            if is_teapot_bend or is_linac_bend:
                if node.getParam("ea1") != 0.0 or node.getParam("ea2") != 0.0:
                    message = f"Found bend ea1 or ea2 != 0.0 ({node.getName()}.)"
                    message += " Nonzero edge angles are not yet supported in envelope tracking."
                    message += " Please set them to zero:"
                    message += "   `node.setParam('ea1', 0.0)`"
                    message += "   `node.setParam('ea2', 0.0)`"
                    raise RuntimeError(message)

    def _getEnvelopeSpaceChargeMatrix(self, envelope: Envelope, length: float, sc: str | None) -> np.ndarray | None:
        if not sc or length <= 0:
            return None
        if sc == "2d":
            return envelope.sc_matrix_2d(length)
        if sc == "3d":
            return envelope.sc_matrix_3d(length)
        raise ValueError(f"Invalid envelope space charge option `{sc}`")

    def setEnvelopeSpaceCharge(self, sc: str | None) -> None:
        self.__envelopeSpaceCharge = sc

    def trackEnvelope(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        history: bool = False,
    ) -> None | dict[str, list]:
        """
        Track envelope through lattice.
        """
        if history:
            return self.trackEnvelopeHistory(
                envelope,
                index_start=index_start,
                index_stop=index_stop,
                sc=sc
            )

        self._prepareEnvelopeTracking()
        self.setEnvelopeSpaceCharge(sc)
        sync_part = envelope.sync_part

        for node in self._getNodesInRange(index_start, index_stop):
            for child_node in node.getChildNodes(AccNode.ENTRANCE):
                matrix = child_node.getMatrix(sync_part)
                if matrix is not None:
                    envelope.transform(matrix)

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(AccNode.BODY, part_index, place_in_part=AccNode.BEFORE):
                    matrix = child_node.getMatrix(sync_part)
                    if matrix is not None:
                        envelope.transform(matrix)

                matrix_sc = self._getEnvelopeSpaceChargeMatrix(envelope, node.getLength(part_index), sc)
                matrix = node.getMatrix(sync_part, part_index=part_index)
                if matrix is not None:
                    if matrix_sc is not None:
                        matrix = matrix @ matrix_sc
                    envelope.transform(matrix)

                for child_node in node.getChildNodes(AccNode.BODY, part_index, place_in_part=AccNode.AFTER):
                    matrix = child_node.getMatrix(sync_part)
                    if matrix is not None:
                        envelope.transform(matrix)

            for child_node in node.getChildNodes(AccNode.EXIT):
                matrix = child_node.getMatrix(sync_part)
                if matrix is not None:
                    envelope.transform(matrix)

    def trackEnvelopeHistory(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
    ) -> dict[str, list]:
        """
        Track envelope and return parameters vs. position in lattice.
        """
        self._prepareEnvelopeTracking()
        self.setEnvelopeSpaceCharge(sc)
        sync_part = envelope.sync_part

        history_keys = [
            "s",
            "kin_energy",
            "gamma",
            "beta",
            "mean",
            "cov",
            "rms_x",
            "rms_y",
            "rms_z",
            "eps_x",
            "eps_y",
        ]
        history = {key: [] for key in history_keys}

        def observe(envelope: Envelope) -> dict:
            parameters = {}
            parameters["gamma"] = envelope.gamma
            parameters["beta"] = envelope.beta
            parameters["kin_energy"] = envelope.kin_energy
            parameters["mean"] = envelope.centroid.copy()
            parameters["cov"] = envelope.cov_matrix.copy()
            parameters["rms_x"] = np.sqrt(parameters["cov"][0, 0])
            parameters["rms_y"] = np.sqrt(parameters["cov"][2, 2])
            parameters["rms_z"] = np.sqrt(parameters["cov"][4, 4])
            return parameters

        def update_history(envelope: Envelope, position: float) -> None:
            history["s"].append(position)
            parameters = observe(envelope)
            for key in parameters:
                history[key].append(parameters[key])

        path_length = 0.0
        update_history(envelope, path_length)

        for node in self._getNodesInRange(index_start, index_stop):
            for child_node in node.getChildNodes(AccNode.ENTRANCE):
                matrix = child_node.getMatrix(sync_part)
                if matrix is not None:
                    envelope.transform(matrix)

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(AccNode.BODY, part_index, place_in_part=AccNode.BEFORE):
                    matrix = child_node.getMatrix(sync_part)
                    if matrix is not None:
                        envelope.transform(matrix)

                matrix_sc = self._getEnvelopeSpaceChargeMatrix(envelope, node.getLength(part_index), sc)
                matrix = node.getMatrix(sync_part, part_index=part_index)
                if matrix is not None:
                    if matrix_sc is not None:
                        matrix = matrix @ matrix_sc
                    envelope.transform(matrix)

                path_length += node.getLength(part_index)
                update_history(envelope, path_length)

                for child_node in node.getChildNodes(AccNode.BODY, part_index, place_in_part=AccNode.AFTER):
                    matrix = child_node.getMatrix(sync_part)
                    if matrix is not None:
                        envelope.transform(matrix)

            for child_node in node.getChildNodes(AccNode.EXIT):
                matrix = child_node.getMatrix(sync_part)
                if matrix is not None:
                    envelope.transform(matrix)
        return history

    def precomputeEnvelopeMatrices(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
    ) -> list:
        """
        Pre-compute transfer matrices for each node.

        For each node, store tuple (node, matrix). Space charge kicks are
        stored as ("sc", length).
        """
        self._prepareEnvelopeTracking()
        sync_part = envelope.sync_part

        self.__envelopeElements = []
        self.__envelopeOneTurnMatrix = None
        self.__envelopeSpaceCharge = sc

        for node in self._getNodesInRange(index_start, index_stop):
            for child_node in node.getChildNodes(AccNode.ENTRANCE):
                matrix = child_node.getMatrix(sync_part)
                if matrix is not None:
                    self.__envelopeElements.append((child_node, matrix))

            for part_index in range(node.getnParts()):
                for child_node in node.getChildNodes(AccNode.BODY, part_index, place_in_part=AccNode.BEFORE):
                    matrix = child_node.getMatrix(sync_part)
                    if matrix is not None:
                        self.__envelopeElements.append((child_node, matrix))

                if sc:
                    length = node.getLength(part_index)
                    if length > 0:
                        self.__envelopeElements.append(("sc", length))

                matrix = node.getMatrix(sync_part, part_index=part_index)
                if matrix is not None:
                    self.__envelopeElements.append((node, matrix))

                for child_node in node.getChildNodes(AccNode.BODY, part_index, place_in_part=AccNode.AFTER):
                    matrix = child_node.getMatrix(sync_part)
                    if matrix is not None:
                        self.__envelopeElements.append((child_node, matrix))

            for child_node in node.getChildNodes(AccNode.EXIT):
                matrix = child_node.getMatrix(sync_part)
                if matrix is not None:
                    self.__envelopeElements.append((child_node, matrix))

        return self.__envelopeElements

    def trackEnvelopeRing(self, envelope: Envelope, sc: str | None = None) -> None:
        """
        Track using pre-computed transfer matrices.

        The method assumes that all nodes are static and that there is no
        change in the synchronous particle energy. In this case the matrices
        can be computed once and reused on each turn. If there is no space charge,
        we track using the one-turn matrix.
        """
        if not self.__envelopeElements or self.__envelopeSpaceCharge != sc:
            self.precomputeEnvelopeMatrices(envelope, sc=sc)

        if not sc:
            if self.__envelopeOneTurnMatrix is None:
                self.__envelopeOneTurnMatrix = np.identity(7)
                for node, matrix in self.__envelopeElements:
                    self.__envelopeOneTurnMatrix = matrix @ self.__envelopeOneTurnMatrix
            envelope.transform(self.__envelopeOneTurnMatrix)
            return

        for element in self.__envelopeElements:
            if element[0] == "sc":
                length = element[1]
                matrix = self._getEnvelopeSpaceChargeMatrix(envelope, length, sc)
                envelope.transform(matrix)
            else:
                node, matrix = element
                envelope.transform(matrix)

    def getEnvelopeTransferMatrix(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
    ) -> np.ndarray:
        """
        Return total transfer matrix, including linear space charge when requested.
        """
        elements = self.precomputeEnvelopeMatrices(envelope, index_start, index_stop, sc=sc)

        total_matrix = np.identity(7)
        for element in elements:
            if element[0] == "sc":
                length = element[1]
                matrix = self._getEnvelopeSpaceChargeMatrix(envelope, length, sc)
            else:
                node, matrix = element
            envelope.transform(matrix)
            total_matrix = matrix @ total_matrix
        return total_matrix
