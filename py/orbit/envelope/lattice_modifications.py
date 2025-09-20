import collections

from ..lattice import AccActionsContainer
from ..lattice import AccLattice
from ..lattice import AccNode
from ..lattice import AccNodeBunchTracker

from .nodes import EnvelopeTrackerNode
from .nodes import KVEnvelopeTrackerNode
from .nodes import DanilovEnvelopeTrackerNode


class Parent:
    def __init__(self, node: AccNode, part_index: int, position: float, path_length: float) -> None:
        self.node = node
        self.name = self.node.getName()
        self.part_index = part_index
        self.position = position
        self.path_length = path_length


def set_max_path_length(lattice: AccLattice, length: float) -> AccLattice:
    if length:
        for node in lattice.getNodes():
            if node.getLength() > length:
                node.setnParts(1 + int(node.getLength() / length))
    return lattice


def add_envelope_tracker_nodes(
    lattice: AccLattice,
    path_length_max: float,
    path_length_min: float,
    constructor: EnvelopeTrackerNode,
    constructor_kwargs: dict,
) -> list[EnvelopeTrackerNode]:

    nodes = lattice.getNodes()
    if not nodes:
        return

    lattice = set_max_path_length(lattice, path_length_max)

    parents = []
    length_total = 0.0
    length_total = running_path = rest_length = 0.0
    for node in nodes:
        for part_index in range(node.getnParts()):
            part_length = node.getLength(part_index)
            parent = Parent(node, part_index, position=length_total, path_length=running_path)
            if running_path > path_length_min:
                parents.append(parent)
                running_path = 0.0
            running_path += part_length
            length_total += part_length

    if len(parents) > 0:
        rest_length = length_total - parents[-1].position
    else:
        rest_length = length_total

    parents.insert(0, Parent(node=nodes[0], part_index=0, position=0.0, path_length=rest_length))

    tracker_nodes = []
    for i in range(len(parents) - 1):
        parent = parents[i]
        parent_new = parents[i + 1]

        tracker_node_name = "{}:{}:".format(parent.name, parent.part_index)
        tracker_node = constructor(
            name=tracker_node_name,
            kick_length=parent_new.path_length,
            **constructor_kwargs,
        )
        parent.node.addChildNode(
            tracker_node, parent.node.BODY, parent.part_index, parent.node.BEFORE
        )
        tracker_nodes.append(tracker_node)

    parent = parents[-1]
    tracker_node = constructor(
        name="{}:{}:".format(parent.node.getName(), parent.part_index),
        kick_length=rest_length,
        **constructor_kwargs,
    )
    tracker_nodes.append(tracker_node)
    parent.node.addChildNode(tracker_node, parent.node.BODY, parent.part_index, parent.node.BEFORE)

    return tracker_nodes


def add_kv_envelope_tracker_nodes(
    lattice: AccLattice, 
    path_length_max: float = None, 
    path_length_min: float = 1.00e-06, 
    **kwargs
) -> None:
    tracker_nodes = add_envelope_tracker_nodes(
        lattice=lattice,
        path_length_max=path_length_max,
        path_length_min=path_length_min,
        constructor=KVEnvelopeTrackerNode,
        constructor_kwargs=kwargs,
    )
    for tracker_node in tracker_nodes:
        name = "".join([tracker_node.getName(), ":", "kv_envelope_tracker"])
        tracker_node.setName(name)
    lattice.initialize()
    return tracker_nodes


def add_danilov_envelope_tracker_nodes(
    lattice: AccLattice, 
    path_length_max: float = None, 
    path_length_min: float = 1.00e-06, 
    **kwargs
) -> None:
    tracker_nodes = add_envelope_tracker_nodes(
        lattice=lattice,
        path_length_max=path_length_max,
        path_length_min=path_length_min,
        constructor=DanilovEnvelopeTrackerNode,
        constructor_kwargs=kwargs,
    )
    for tracker_node in tracker_nodes:
        name = "".join([tracker_node.getName(), ":", "danilov_envelope_tracker"])
        tracker_node.setName(name)
    lattice.initialize()
    return tracker_nodes