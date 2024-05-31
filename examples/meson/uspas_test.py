# This script runs through how to build a lattice, build a bunch, and how to track it.

import numpy as np
from matplotlib import pyplot as plt

# Import Bunch to be able to build a bunch.
from orbit.core.bunch import Bunch

# Import LinacAccLattice to build a lattice and Sequence to build a sequence that we can then import into the lattice.
from orbit.py_linac.lattice.LinacAccLatticeLib import LinacAccLattice, Sequence

# Import Quad and Drift components to build the lattice.
from orbit.py_linac.lattice.LinacAccNodes import Quad, Drift

# Import AccActionsContainer, a method to add functionality throughout the accelerator.
from orbit.lattice import AccActionsContainer

# Field strength and length of the quadrupoles
field_str = 1.8
quad_len = 0.2

# Lengths of the interior and exterior drifts.
drift_end_len = 0.4
drift_len = 0.8

# List to contain the drift and quadrupole nodes.
list_of_nodes = []

# What follows is defining each node of the lattice.
# The first drift node (needs a specified length)
D1 = Drift("Drift1")
D1.setLength(drift_end_len)
list_of_nodes.append(D1)

# The first quadrupole node (needs a specified field gradient strength and length)
Q1 = Quad("Quad1")
Q1.setLength(quad_len)
Q1.setField(field_str)
list_of_nodes.append(Q1)

# The second drift node
D2 = Drift("Drift2")
D2.setLength(drift_len)
list_of_nodes.append(D2)

# The second quadrupole node
Q2 = Quad("Quad2")
Q2.setLength(quad_len)
Q2.setField(-field_str)
list_of_nodes.append(Q2)

# The third drift node
D3 = Drift("Drift3")
D3.setLength(drift_len)
list_of_nodes.append(D3)

# The third quadrupole node
Q3 = Quad("Quad3")
Q3.setLength(quad_len)
Q3.setField(field_str)
list_of_nodes.append(Q3)

# The fourth drift node
D4 = Drift("Drift4")
D4.setLength(drift_len)
list_of_nodes.append(D4)

# The fourth quadrupole node
Q4 = Quad("Quad4")
Q4.setLength(quad_len)
Q4.setField(-field_str)
list_of_nodes.append(Q4)

# The fifth drift node
D5 = Drift("Drift5")
D5.setLength(drift_end_len)
list_of_nodes.append(D5)

# Define the sequence and add the list of nodes to the sequence.
fodo = Sequence('FODO')
fodo.setNodes(list_of_nodes)

# Define the lattice, add the list of nodes to the lattice, and initialize the lattice.
my_lattice = LinacAccLattice('My Lattice')
my_lattice.setNodes(list_of_nodes)
my_lattice.initialize()
print("Total length=", my_lattice.getLength())

# Define the bunch and it's kinetic energy.
bunch = Bunch()
bunch.getSyncParticle().kinEnergy(0.0025)  # GeV

# Add macroparticles to the bunch.
bunch.addParticle(0.00, 0.002, -0.001, 0.001, 0.0, 0.0)
bunch.addParticle(0.001, 0.000, 0.001, -0.001, 0.0, 0.0)
bunch.addParticle(-0.001, -0.002, 0.001, -0.001, 0.0, 0.0)
for i in range(100000):
    bunch.addParticle(-0.001, -0.002, 0.001, -0.001, 0.0, 0.0)

# Setup RF cavities (if we had any).
my_lattice.trackDesignBunch(bunch)

# Prepare variables for use in our track action. We have an array of the position along the linac and an array of x
# positions.
pos_array = []
x_array = []

# Create a parameter dictionary to pass to our track action.
my_params = {"old_pos": -1.0}


# Action function to track the horizontal position of the particle through the lattice. This will be used at the
# entrance of each node.
def action_entrance(paramsDict):
    # Grab values from the parameter dictionary.
    node = paramsDict["node"]
    bunch = paramsDict["bunch"]
    pos = paramsDict["path_length"]

    # Don't continue if current node has the same position as the previous node.
    if pos <= paramsDict["old_pos"]:
        return

    # Replace the old_pos value with the current nodes position, then add the position to out position array.
    paramsDict["old_pos"] = pos
    pos_array.append(pos)

    # Loop through the particles and add their x positions to a list, then add that list to our x_array.
    nParts = bunch.getSizeGlobal()
    x_temp = []
    for n in range(nParts):
        x = bunch.x(n) * 1000  # [mm]
        x_temp.append(x)
    x_array.append(x_temp)


# Use the same function for an action that will take place at the exit of each node.
def action_exit(paramsDict):
    action_entrance(paramsDict)


# Define your action container and add your actions. One action is added to the entrance of each node, the other to the
# exit of each node.
my_container = AccActionsContainer("Test Bunch Tracking")
my_container.addAction(action_entrance, AccActionsContainer.ENTRANCE)
my_container.addAction(action_exit, AccActionsContainer.EXIT)

# Track the bunch through the lattice. We are also passing our parameter dictionary and action container.
my_lattice.trackBunch(bunch, paramsDict=my_params, actionContainer=my_container)

# Convert x_array into a numpy array. Then loop through each particle to plot their positions through the lattice.
x_array = np.array(x_array)
num_of_parts = bunch.getSizeGlobal()
for n in range(3):
    plt.plot(pos_array, x_array[:, n], label='Particle ' + str(n+1))
plt.xlabel('Lattice position [m]')
plt.ylabel('Horizontal Position [mm]')
plt.legend()
plt.show()