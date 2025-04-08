import collections

from ..lattice import AccActionsContainer
from ..lattice import AccLattice
from ..lattice import AccNode
from ..lattice import AccNodeBunchTracker

from .danilov_envelope_solver_nodes import DanilovEnvelopeSolverNode20
from .danilov_envelope_solver_nodes import DanilovEnvelopeSolverNode22


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


def add_danilov_envolope_solver_nodes(
    lattice: AccLattice, 
    path_length_max: float,
    path_length_min: float, 
    solver_node_constructor: DanilovEnvelopeSolverNode20 | DanilovEnvelopeSolverNode22, 
    solver_node_constructor_kwargs: dict,
) -> list[DanilovEnvelopeSolverNode20 | DanilovEnvelopeSolverNode22]:
    
    nodes = lattice.getNodes()
    if not nodes:
        return
    
    lattice = set_max_path_length(lattice, path_length_max)
    
    parents = []
    length_total = running_path = rest_length = 0.0
    for node in nodes:
        for part_index in range(node.getnParts()):
            part_length = node.getLength(part_index)
            if part_length > 1.0:
                message  = "Warning! Node {} has length {} > 1 m. ".format(node.getName(), part_length)
                message += "Space charge algorithm may be innacurate!"
                print(message)

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
    
    solver_nodes = []
    for i in range(len(parents) - 1):
        parent = parents[i]
        parent_new = parents[i + 1]

        solver_node_name = "{}:{}:".format(parent.name, parent.part_index)
        solver_node = solver_node_constructor(
            name=solver_node_name,
            kick_length=parent_new.path_length,
            **solver_node_constructor_kwargs
        )
        parent.node.addChildNode(solver_node, parent.node.BODY, parent.part_index, parent.node.BEFORE)
        solver_nodes.append(solver_node)
        
    parent = parents[-1]
    solver_node = solver_node_constructor(
        name="{}:{}:".format(parent.node.getName(), parent.part_index),
        kick_length=rest_length,
        **solver_node_constructor_kwargs
    )
    solver_nodes.append(solver_node)
    parent.node.addChildNode(solver_node, parent.node.BODY, parent.part_index, parent.node.BEFORE)

    return solver_nodes


def add_danilov_envelope_solver_nodes_20(
    lattice: AccLattice, 
    path_length_max: float = None,
    path_length_min: float = 1.00e-06,   
    **kwargs  
) -> None:
    solver_nodes = add_danilov_envolope_solver_nodes(
        lattice=lattice, 
        path_length_max=path_length_max,
        path_length_min=path_length_min, 
        solver_node_constructor=DanilovEnvelopeSolverNode20, 
        solver_node_constructor_kwargs=kwargs
    )
    for solver_node in solver_nodes:
        name = "".join([solver_node.getName(), ":", "danilov_env_solver_20"])
        solver_node.setName(name)
    lattice.initialize()
    return solver_nodes


def add_danilov_envelope_solver_nodes_22(
    lattice: AccLattice, 
    path_length_max: float = None,
    path_length_min: float = 1.00e-06,   
    **kwargs  
) -> None:
    solver_nodes = add_danilov_envolope_solver_nodes(
        lattice=lattice, 
        path_length_max=path_length_max,
        path_length_min=path_length_min, 
        solver_node_constructor=DanilovEnvelopeSolverNode22, 
        solver_node_constructor_kwargs=kwargs
    )
    for solver_node in solver_nodes:
        name = "".join([solver_node.getName(), ":", "danilov_env_solver_22"])
        solver_node.setName(name)
    lattice.initialize()
    return solver_nodes