import os

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
        self._envelope_matrices = {}
        self._envelope_total_matrix = None
        self._envelope_cache_key = None
        self._envelope_cache_time = None

    def initialize(self):
        """
        Method. Initializes the lattice and child node structures.
        """
        self._envelope_cache_key = None
        self._envelope_total_matrix = None
        self._envelope_cache_time = None
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

    def trackEnvelope(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        paramsDict: dict | None = None,
        actionContainer: AccActionsContainer | None = None,
        fit: bool = True,
        static: bool = False,
        fallback: bool = True,
    ) -> None:
        """Track envelope through the lattice.

        Args:
            envelope: Envelope to track.
            index_start: Index of first node in sublattice.
            index_stop: Index of last node in sublattice.
            sc: Whether to include space charge kicks.
            paramsDict: Shared action context. Populated with envelope, lattice,
                actions, node, parentNode, and path_length (default 0). An existing
                path_length is continued, as in bunch tracking.
            actionContainer: AccActionsContainer with user actions. ENTRANCE runs
                before entrance children; BODY runs before each part's space
                charge and transport; EXIT runs after exit children. Actions also
                run on child nodes, in registration order. Return values are ignored.
            fit: Whether to fit all transfer matrices. With False, use analytic
                matrices, except bends with nonzero ea1 or ea2 (including their
                children), which always use fitted matrices. This is set True
                by default for more accurate simulations with nonlinear elements,
                but it should be turned off to increase speed in linear lattices.
            fallback: With fit=False, fit nodes whose getMatrix raises
                NotImplementedError. A None matrix is an intentional no-op.
            static: Whether to cache transfer matrices during tracking and
                reuse them for the same incoming synchronous-particle state.
                Assumes the lattice has no time-dependent or accelerating nodes.
                Synchronous-particle time still advances on each traversal.
                Changes to node parameters require reinitializing the lattice
                to invalidate the cache.
        """
        self._trackEnvelope(
            envelope,
            index_start,
            index_stop,
            sc,
            paramsDict,
            actionContainer,
            fit,
            static,
            fallback,
        )

    def _trackEnvelope(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int | None = None,
        sc: str | None = None,
        paramsDict: dict | None = None,
        actionContainer: AccActionsContainer | None = None,
        fit: bool = True,
        static: bool = False,
        fallback: bool = True,
        calculate_matrix: bool = False,
    ):
        from orbit.py_linac.lattice.LinacAccNodes import Bend
        from orbit.teapot.teapot import BendTEAPOT

        if sc not in (None, False, "2d", "3d"):
            raise ValueError(f"Invalid envelope space charge option `{sc}`")

        params = {} if paramsDict is None else paramsDict
        actions = actionContainer if actionContainer is not None else AccActionsContainer()

        params["envelope"] = envelope
        params["lattice"] = self
        params["actions"] = actions
        params.setdefault("path_length", 0.0)
        if hasattr(self, "getUseRealCharge"):
            params.setdefault("useCharge", self.getUseRealCharge())

        start_position = params["path_length"]
        sync_part = envelope.sync_part
        combine = static and not sc and not any(actions.getActions(p) for p in (0, 1, 2))
        total_matrix = np.eye(7) if calculate_matrix or combine else None

        start, stop, _ = slice(index_start, None if index_stop is None else index_stop + 1).indices(len(self.__children))

        if static:
            cache_key = (
                start,
                stop,
                fit,
                fallback,
                sync_part.mass(),
                sync_part.charge(),
                sync_part.rVector(),
                sync_part.pVector(),
                sync_part.nxVector(),
            )
            if self._envelope_cache_key != cache_key:
                self._envelope_cache_key = cache_key
                self._envelope_matrices = {}
                self._envelope_total_matrix = None
                self._envelope_cache_time = sync_part.time()
            time_offset = sync_part.time() - self._envelope_cache_time
            if combine and self._envelope_total_matrix is not None:
                envelope.transform(self._envelope_total_matrix)
                reference, distance, last_node, last_parent = self._envelope_end_state
                copy_sync_particle(reference.getSyncParticle(), sync_part)
                sync_part.time(sync_part.time() + time_offset)
                params["path_length"] += distance
                if last_node is not None:
                    params["node"], params["parentNode"] = last_node, last_parent
                return self._envelope_total_matrix if calculate_matrix else None

        probe_bunch = None
        generator = None
        fit_stack = [fit]

        def enter(context):
            node = context["node"]
            edge_bend = isinstance(node, (Bend, BendTEAPOT)) and (
                node.getParam("ea1") != 0.0 or node.getParam("ea2") != 0.0
            )
            fit_stack.append(fit_stack[-1] or edge_bend)

        def exit(context):
            fit_stack.pop()

        def apply(matrix):
            nonlocal total_matrix
            if matrix is not None:
                envelope.transform(matrix)
                if total_matrix is not None:
                    total_matrix = matrix @ total_matrix

        def transport(context):
            nonlocal probe_bunch, generator
            node = context["node"]
            part = node.getActivePartIndex()
            length = node.getLength(part)
            if sc and length > 0:
                apply(
                    envelope.sc_matrix_2d(length) if sc == "2d" else envelope.sc_matrix_3d(length)
                )

            key = (node, part)
            if static and key in self._envelope_matrices:
                matrix, reference = self._envelope_matrices[key]
                copy_sync_particle(reference.getSyncParticle(), sync_part)
                sync_part.time(sync_part.time() + time_offset)
            else:
                use_fit = fit_stack[-1]
                if not use_fit:
                    try:
                        matrix = node.getMatrix(sync_part, part_index=part)
                    except NotImplementedError:
                        if not fallback:
                            raise
                        use_fit = True
                if use_fit:
                    if probe_bunch is None:
                        probe_bunch = bunch_from_sync_particle(sync_part)
                        generator = MatrixGenerator()
                    copy_sync_particle(sync_part, probe_bunch.getSyncParticle())
                    matrix = fit_node_transfer_matrix(
                        node,
                        probe_bunch,
                        part_index=part,
                        parent_node=context["parentNode"],
                        matrix_generator=generator,
                        params_dict=context,
                    )
                    copy_sync_particle(probe_bunch.getSyncParticle(), sync_part)
                if static:
                    self._envelope_matrices[key] = (
                        matrix,
                        bunch_from_sync_particle(sync_part),
                    )
            apply(matrix)

        # Append transport after user BODY actions, exactly as trackBunch does.
        internal_actions = (
            (enter, AccNode.ENTRANCE),
            (transport, AccNode.BODY),
            (exit, AccNode.EXIT),
        )
        for action, place in internal_actions:
            actions.addAction(action, place)
        try:
            if start < stop:
                self.trackActions(actions, params, start, stop - 1)
        finally:
            for action, place in internal_actions:
                actions.removeAction(action, place)

        if combine:
            self._envelope_total_matrix = total_matrix
            self._envelope_end_state = (
                bunch_from_sync_particle(sync_part),
                params["path_length"] - start_position,
                params.get("node"),
                params.get("parentNode"),
            )
        return total_matrix

    def getEnvelopeTransferMatrix(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        fit: bool = True,
        fallback: bool = False,
    ) -> np.ndarray:
        return self._trackEnvelope(
            envelope.copy(),
            index_start,
            index_stop,
            sc,
            fit=fit,
            fallback=fallback,
            calculate_matrix=True,
        )
