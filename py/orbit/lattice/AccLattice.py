import os
from typing import Any
from typing import Callable

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle
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

    def _prepareEnvelopeTracking(self, index_start: int, index_stop: int, fit: bool = True) -> None:
        """Check lattice before tracking envelope."""
        if fit:
            return

        from orbit.py_linac.lattice.LinacAccNodes import Bend
        from orbit.teapot.teapot import BendTEAPOT

        for node in self.__children:
            if isinstance(node, BendTEAPOT):
                uses_unsupported_fringe = (
                    node.getParam("ea1") != 0.0 and node.getUsageFringeFieldIN()
                ) or (node.getParam("ea2") != 0.0 and node.getUsageFringeFieldOUT())
                if uses_unsupported_fringe:
                    message = f"Found an enabled fringe field with a nonzero edge angle ({node.getName()})."
                    message += (
                        " Analytic envelope tracking supports the wedge transformations only."
                    )
                    message += " Disable the bend fringe field or use `fit=True`."
                    raise RuntimeError(message)
            if isinstance(node, Bend):
                if node.getParam("ea1") != 0.0 or node.getParam("ea2") != 0.0:
                    message = f"Found bend ea1 or ea2 != 0.0 ({node.getName()}.)"
                    message += " Nonzero edge angles are not yet supported in envelope tracking."
                    message += " Please set them to zero:"
                    message += "   `node.setParam('ea1', 0.0)`"
                    message += "   `node.setParam('ea2', 0.0)`"
                    raise RuntimeError(message)

    @staticmethod
    def _getEnvelopeCacheKey(envelope: Envelope, index_start: int, index_stop: int, sc: bool, fit: bool):
        sync_part = envelope.sync_part
        sync_state = (
            sync_part.kinEnergy(),
            sync_part.time(),
            sync_part.mass(),
            sync_part.charge(),
        )
        return index_start, index_stop, sc, fit, sync_state

    def _createEnvelopeFitState(self, sync_part) -> dict[str, Any]:
        """Return items needed to compute best-fit transfer matrix."""
        bunch = bunch_from_sync_particle(sync_part)
        lost_bunch = Bunch()
        bunch.copyEmptyBunchTo(lost_bunch)
        params_dict = {}
        if hasattr(self, "getUseRealCharge"):
            params_dict["useCharge"] = self.getUseRealCharge()
        return {
            "bunch": bunch,
            "lost_bunch": lost_bunch,
            "matrix_generator": MatrixGenerator(),
            "params_dict": params_dict,
        }

    def _getEnvelopeNodeMatrix(
        self,
        node: AccNode,
        sync_part: SyncParticle,
        part_index: int | None = None,
        parent_node: AccNode = None,
        fit_state: dict | None = None,
    ) -> np.ndarray | None:
        """Return 7 x 7 transfer matrix from node and synchronous particle.

        If `fit_state` is provided, the matrix will be fit to input/output
        coordinates of a particle launched near the origin. Otherwise the
        analytic matrix will be used (if implemented).
        """
        if fit_state is None:
            if part_index is None:
                return node.getMatrix(sync_part)
            return node.getMatrix(sync_part, part_index=part_index)

        matrix = fit_node_transfer_matrix(
            node,
            fit_state["bunch"],
            part_index=0 if part_index is None else part_index,
            parent_node=parent_node,
            matrix_generator=fit_state["matrix_generator"],
            lost_bunch=fit_state["lost_bunch"],
            params_dict=fit_state["params_dict"],
        )
        copy_sync_particle(fit_state["bunch"].getSyncParticle(), sync_part)
        return matrix

    def _getEnvelopeSpaceChargeMatrix(
        self, envelope: Envelope, length: float, sc: str | None
    ) -> np.ndarray | None:
        """Return transfer matrix for linear space charge kick."""
        if not sc or length <= 0:
            return None
        if sc == "2d":
            return envelope.sc_matrix_2d(length)
        if sc == "3d":
            return envelope.sc_matrix_3d(length)
        raise ValueError(f"Invalid envelope space charge option `{sc}`")

    def _iterateEnvelopeElements(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        fit: bool = True,
    ):
        """Yield matrices, space charge kicks, and position updates in tracking order."""
        self._prepareEnvelopeTracking(index_start, index_stop, fit=fit)
        sync_part = envelope.sync_part
        fit_state = self._createEnvelopeFitState(sync_part) if fit else None

        def iter_child_elements(child_nodes, parent_node):
            for child_node in child_nodes:
                matrix = self._getEnvelopeNodeMatrix(
                    child_node,
                    sync_part,
                    parent_node=parent_node,
                    fit_state=fit_state,
                )
                if matrix is not None:
                    yield child_node, matrix

        for node in self._getNodesInRange(index_start, index_stop):
            yield from iter_child_elements(
                node.getChildNodes(AccNode.ENTRANCE),
                node,
            )

            for part_index in range(node.getnParts()):
                yield from iter_child_elements(
                    node.getChildNodes(
                        AccNode.BODY,
                        part_index,
                        place_in_part=AccNode.BEFORE,
                    ),
                    node,
                )

                length = node.getLength(part_index)
                if sc and length > 0:
                    yield "sc", length

                matrix = self._getEnvelopeNodeMatrix(
                    node,
                    sync_part,
                    part_index=part_index,
                    parent_node=self,
                    fit_state=fit_state,
                )
                if matrix is not None:
                    yield node, matrix

                yield "position", length

                yield from iter_child_elements(
                    node.getChildNodes(
                        AccNode.BODY,
                        part_index,
                        place_in_part=AccNode.AFTER,
                    ),
                    node,
                )

            yield from iter_child_elements(
                node.getChildNodes(AccNode.EXIT),
                node,
            )

    def _applyEnvelopeElements(
        self,
        envelope: Envelope,
        elements: list,
        sc: str | None = None,
        update_history: Callable = None,
        calculate_matrix: bool = False,
    ) -> np.ndarray | None:
        """Apply envelope operations and optionally return their combined matrix."""
        total_matrix = np.identity(7) if calculate_matrix else None
        path_length = 0.0

        for element_type, value in elements:
            if element_type == "position":
                path_length += value
                if update_history is not None:
                    update_history(path_length)
                continue

            if element_type == "sc":
                matrix = self._getEnvelopeSpaceChargeMatrix(envelope, value, sc)
            else:
                matrix = value

            envelope.transform(matrix)
            if calculate_matrix:
                total_matrix = matrix @ total_matrix

        return total_matrix

    def _precomputeEnvelopeElements(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        fit: bool = True,
    ) -> list:
        """Precompute envelope elements in lattice."""
        cache_key = self._getEnvelopeCacheKey(envelope, index_start, index_stop, sc, fit)
        self._envelope_elements = list(
            self._iterateEnvelopeElements(
                envelope,
                index_start=index_start,
                index_stop=index_stop,
                sc=sc,
                fit=fit,
            )
        )
        self._envelope_cache_key = cache_key
        self._envelope_total_matrix = None
        return self._envelope_elements

    def _getStaticEnvelopeElements(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        fit: bool = True,
    ) -> list:
        cache_key = self._getEnvelopeCacheKey(envelope, index_start, index_stop, sc, fit)
        if self._envelope_cache_key != cache_key:
            self._precomputeEnvelopeElements(
                envelope,
                index_start=index_start,
                index_stop=index_stop,
                sc=sc,
                fit=fit,
            )
        return self._envelope_elements

    def trackEnvelope(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        history: bool = False,
        fit: bool = True,
        static: bool = False,
    ) -> None | dict[str, list]:
        """Track envelope through the lattice.

        Args:
            envelope: Envelope to track.
            index_start: Index of first node in sublattice.
            index_stop: Index of last node in sublattice.
            sc: Whether to include space charge kicks.
            history: Whether to return beam parameters vs. position in lattice.
            fit: Whether to use best-fit transfer matrices or analytic
                transfer matrices.
            static: Whether to pre-compute transfer matrices before
                tracking. This works for as long as there are no time-dependent
                nodes in the lattice.
        """
        if history:
            return self._trackEnvelopeHistory(
                envelope,
                index_start=index_start,
                index_stop=index_stop,
                sc=sc,
                fit=fit,
                static=static,
            )
        if static:
            self._trackEnvelopeStatic(
                envelope,
                index_start=index_start,
                index_stop=index_stop,
                sc=sc,
                fit=fit,
            )
            return None

        elements = self._iterateEnvelopeElements(
            envelope,
            index_start=index_start,
            index_stop=index_stop,
            sc=sc,
            fit=fit,
        )
        self._applyEnvelopeElements(envelope, elements, sc=sc)

    def _trackEnvelopeHistory(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        fit: bool = True,
        static: bool = False,
    ) -> dict[str, list]:
        """Track envelope and return parameters vs. position in lattice."""
        keys = [
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
            "eps_x_n",
            "eps_y",
            "eps_y_n",
            "eps_z",
            "eps_z_n",
            "eps_1",
            "eps_2",
            "eps_3",
        ]
        history = {key: [] for key in keys}

        def observe(envelope: Envelope) -> dict:
            cov_matrix = envelope.cov_matrix.copy()
            poisson_matrix = np.zeros_like(cov_matrix)
            for i in range(0, 6, 2):
                poisson_matrix[i, i + 1] = 1.0
                poisson_matrix[i + 1, i] = -1.0

            eigvals = np.linalg.eigvals(cov_matrix @ poisson_matrix)
            eigenemittances = np.imag(eigvals)
            eigenemittances = eigenemittances[eigenemittances > 0]

            gamma = envelope.gamma
            beta = envelope.beta
            emittances = [
                np.sqrt(np.linalg.det(cov_matrix[i : i + 2, i : i + 2])) for i in (0, 2, 4)
            ]
            return {
                "gamma": gamma,
                "beta": beta,
                "kin_energy": envelope.kin_energy,
                "mean": envelope.centroid.copy(),
                "cov": cov_matrix,
                "rms_x": np.sqrt(cov_matrix[0, 0]),
                "rms_y": np.sqrt(cov_matrix[2, 2]),
                "rms_z": np.sqrt(cov_matrix[4, 4]),
                "eps_x": emittances[0],
                "eps_y": emittances[1],
                "eps_z": emittances[2],
                "eps_x_n": emittances[0] * gamma * beta,
                "eps_y_n": emittances[1] * gamma * beta,
                "eps_z_n": emittances[2] / beta,
                "eps_1": eigenemittances[0],
                "eps_2": eigenemittances[1],
                "eps_3": eigenemittances[2],
            }

        def update_history(position: float) -> None:
            history["s"].append(position)
            parameters = observe(envelope)
            for key, value in parameters.items():
                history[key].append(value)

        update_history(0.0)
        if static:
            elements = self._getStaticEnvelopeElements(
                envelope,
                index_start=index_start,
                index_stop=index_stop,
                sc=sc,
                fit=fit,
            )
        else:
            elements = self._iterateEnvelopeElements(
                envelope,
                index_start=index_start,
                index_stop=index_stop,
                sc=sc,
                fit=fit,
            )
        self._applyEnvelopeElements(
            envelope,
            elements,
            sc=sc,
            update_history=update_history,
        )
        return history

    def _trackEnvelopeStatic(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        fit: bool = True,
    ) -> None:
        """
        Track using pre-computed transfer matrices.

        The method assumes that all nodes are static and that there is no
        change in the synchronous particle energy. In this case the matrices
        can be computed once and reused on each turn. If there is no space charge,
        we track using the one-turn matrix.
        """
        elements = self._getStaticEnvelopeElements(
            envelope,
            index_start=index_start,
            index_stop=index_stop,
            sc=sc,
            fit=fit,
        )

        if not sc:
            if self._envelope_total_matrix is None:
                self._envelope_total_matrix = np.identity(7)
                for element_type, matrix in elements:
                    if element_type != "position":
                        self._envelope_total_matrix = matrix @ self._envelope_total_matrix
            envelope.transform(self._envelope_total_matrix)
            return

        self._applyEnvelopeElements(envelope, elements, sc=sc)

    def getEnvelopeTransferMatrix(
        self,
        envelope: Envelope,
        index_start: int = 0,
        index_stop: int = None,
        sc: str | None = None,
        fit: bool = True,
    ) -> np.ndarray:
        """Return total transfer matrix, including linear space charge when requested."""
        envelope_out = envelope.copy()
        elements = self._precomputeEnvelopeElements(
            envelope_out,
            index_start=index_start,
            index_stop=index_stop,
            sc=sc,
            fit=fit,
        )
        return self._applyEnvelopeElements(
            envelope_out,
            elements,
            sc=sc,
            calculate_matrix=True,
        )
