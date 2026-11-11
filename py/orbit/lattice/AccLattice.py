import os
from copy import deepcopy
from typing import Any
from typing import Callable

import numpy as np

from orbit.core.teapot_base import MatrixGenerator
from orbit.envelope import Envelope
from orbit.envelope.matrix_fitting import bunch_from_sync_particle
from orbit.envelope.matrix_fitting import copy_sync_particle
from orbit.envelope.matrix_fitting import fit_node_transfer_matrix

from orbit.utils import orbitFinalize
from orbit.utils import NamedObject
from orbit.utils import TypedObject

from .AccActionsContainer import AccActionsContainer
from .AccNode import AccNode


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
        self._envelope_elements = []
        self._envelope_total_matrix = None
        self._envelope_cache_key = None

    def initialize(self):
        """
        Method. Initializes the lattice and child node structures.
        """
        self._envelope_cache_key = None
        self._envelope_total_matrix = None
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

    def _getNodesInRange(self, index_start: int = 0, index_stop: int = None) -> list[AccNode]:
        if index_stop is None:
            index_stop = len(self.__children) - 1
        return self.__children[index_start : index_stop + 1]

    def _iterateEnvelopeElements(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        fit: bool = True,
        fallback: bool = False,
    ):
        """Yield matrices, space charge kicks, and position updates in tracking order."""
        from orbit.py_linac.lattice.LinacAccNodes import Bend
        from orbit.teapot.teapot import BendTEAPOT

        # Work on a private reference particle so precomputation cannot advance
        # the envelope. Sync events replay its state alongside the matrices.
        bunch = bunch_from_sync_particle(envelope.sync_part)
        sync_part = bunch.getSyncParticle()
        matrix_generator = None
        params = {}
        if hasattr(self, "getUseRealCharge"):
            params["useCharge"] = self.getUseRealCharge()

        def node_elements(node, parent_node, part_index=None, force_fit=False):
            nonlocal matrix_generator
            use_fit = fit or force_fit
            if not use_fit:
                try:
                    if part_index is None:
                        matrix = node.getMatrix(sync_part)
                    else:
                        matrix = node.getMatrix(sync_part, part_index=part_index)
                except NotImplementedError:
                    if not fallback:
                        raise
                    use_fit = True
            if use_fit:
                if matrix_generator is None:
                    matrix_generator = MatrixGenerator()
                matrix = fit_node_transfer_matrix(
                    node, bunch,
                    part_index=0 if part_index is None else part_index,
                    parent_node=parent_node,
                    matrix_generator=matrix_generator,
                    params_dict=params,
                )
            yield "sync", bunch_from_sync_particle(sync_part)
            if matrix is not None:
                yield "matrix", matrix

        def child_elements(child_nodes, parent_node, force_fit):
            for child_node in child_nodes:
                yield from node_elements(child_node, parent_node, force_fit=force_fit)

        for node in self._getNodesInRange(index_start, index_stop):
            force_fit = isinstance(node, (Bend, BendTEAPOT)) and (
                node.getParam("ea1") != 0.0 or node.getParam("ea2") != 0.0
            )
            yield from child_elements(
                node.getChildNodes(AccNode.ENTRANCE),
                node,
                force_fit,
            )

            for part_index in range(node.getnParts()):
                yield from child_elements(
                    node.getChildNodes(
                        AccNode.BODY,
                        part_index,
                        place_in_part=AccNode.BEFORE,
                    ),
                    node,
                    force_fit,
                )

                length = node.getLength(part_index)
                if sc and length > 0:
                    yield "sc", length

                yield from node_elements(node, self, part_index, force_fit)

                yield "position", length

                yield from child_elements(
                    node.getChildNodes(
                        AccNode.BODY,
                        part_index,
                        place_in_part=AccNode.AFTER,
                    ),
                    node,
                    force_fit,
                )

            yield from child_elements(
                node.getChildNodes(AccNode.EXIT),
                node,
                force_fit,
            )

    def _applyEnvelopeElements(
        self,
        envelope: Envelope,
        elements,
        sc: str | None = None,
        observe: Callable[[Envelope, float], None] | None = None,
        calculate_matrix: bool = False,
    ) -> np.ndarray | None:
        """Apply envelope operations and optionally return their combined matrix."""
        total_matrix = np.identity(7) if calculate_matrix else None
        path_length = 0.0

        for element_type, value in elements:
            if element_type == "position":
                path_length += value
                if observe is not None:
                    observe(envelope, path_length)
                continue

            if element_type == "sync":
                copy_sync_particle(value.getSyncParticle(), envelope.sync_part)
                continue
            if element_type == "sc":
                if sc == "2d":
                    matrix = envelope.sc_matrix_2d(value)
                elif sc == "3d":
                    matrix = envelope.sc_matrix_3d(value)
                else:
                    raise ValueError(f"Invalid envelope space charge option `{sc}`")
            else:
                matrix = value

            envelope.transform(matrix)
            if calculate_matrix:
                total_matrix = matrix @ total_matrix

        return total_matrix

    def trackEnvelope(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        observe: Callable[[Envelope, float], dict[str, Any]] | None = None,
        fit: bool = True,
        static: bool = False,
        fallback: bool = False,
    ) -> None | dict[str, list]:
        """Track envelope through the lattice.

        Args:
            envelope: Envelope to track.
            index_start: Index of first node in sublattice.
            index_stop: Index of last node in sublattice.
            sc: Whether to include space charge kicks.
            observe: Optional callback observe(envelope, s) returning a dictionary
                with the same keys at every sample. Called at s=0 before tracking
                and after each body part, before its AFTER/EXIT children, matching
                the tracking order. Positions are relative to the selected range.
                Values are copied and collected into a dictionary of lists;
                include "s" in the returned dictionary to record positions.
                With no observer, return None. Observers should not mutate the
                envelope. Use orbit.envelope.EnvelopeAnalysis for beam quantities.
            fit: Whether to fit all transfer matrices. With False, use analytic
                matrices, except bends with nonzero ea1 or ea2 (including their
                children), which always use fitted matrices.
            fallback: With fit=False, fit nodes whose getMatrix raises
                NotImplementedError. A None matrix is an intentional no-op.
            static: Whether to pre-compute transfer matrices before
                tracking and reuse them for the same synchronous-particle state.
                Requires a static lattice; changes to node parameters require
                reinitializing the lattice to invalidate the cache.
        """
        if sc not in (None, False, "2d", "3d"):
            raise ValueError(f"Invalid envelope space charge option `{sc}`")

        elements = self._iterateEnvelopeElements(
            envelope, index_start, index_stop, sc, fit, fallback,
        )
        if static:
            sync = envelope.sync_part
            cache_key = (
                index_start, index_stop, sc, fit, fallback,
                sync.mass(), sync.charge(), sync.time(),
                sync.rVector(), sync.pVector(), sync.nxVector(),
            )
            if self._envelope_cache_key != cache_key:
                self._envelope_elements = list(elements)
                self._envelope_cache_key = cache_key
                self._envelope_total_matrix = None
            elements = self._envelope_elements
            if not sc and observe is None:
                if self._envelope_total_matrix is None:
                    total = np.identity(7)
                    for kind, value in elements:
                        if kind == "matrix":
                            total = value @ total
                    self._envelope_total_matrix = total
                envelope.transform(self._envelope_total_matrix)
                for kind, value in reversed(elements):
                    if kind == "sync":
                        copy_sync_particle(value.getSyncParticle(), sync)
                        break
                return None

        records = {} if observe is not None else None
        record_keys = None

        def record(envelope: Envelope, s: float) -> None:
            nonlocal record_keys
            parameters = observe(envelope, s)
            if record_keys is None:
                record_keys = set(parameters)
            elif parameters.keys() != record_keys:
                raise ValueError("Envelope observer must return the same keys at every sample")
            for key, value in parameters.items():
                records.setdefault(key, []).append(deepcopy(value))

        if observe is not None:
            record(envelope, 0.0)
        self._applyEnvelopeElements(
            envelope, elements, sc=sc,
            observe=record if observe is not None else None,
        )
        return records

    def getEnvelopeTransferMatrix(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        fit: bool = True,
        fallback: bool = False,
    ) -> np.ndarray:
        envelope_out = envelope.copy()
        elements = self._iterateEnvelopeElements(
            envelope_out,
            index_start=index_start,
            index_stop=index_stop,
            sc=sc,
            fit=fit,
            fallback=fallback,
        )
        return self._applyEnvelopeElements(
            envelope_out,
            elements,
            sc=sc,
            calculate_matrix=True,
        )
