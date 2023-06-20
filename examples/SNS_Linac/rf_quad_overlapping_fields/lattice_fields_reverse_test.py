#! /usr/bin/env python

"""
This script compares the gradients, their derivatives, and electric field fo RF gaps
for a lattice and the reversed lattice. The longitunal position for the reversed
lattice is calculated from the end to the start, so plots should be the same.
The lattice can be 4 different types
0. To transformation - hard edge quads and thin RF gaps
1. Replace thin RF-gap node with distributed RF Ez(z) field on the axis
2. Replace hard-edge quads with soft-edge ones, and the quads' fields could overlap
3. Replace RF-gap thin nodes and hard-edge quads with nodes with distributed fields

Script uses Matplotlib and Numpy modules.
"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_to_AxisField_Nodes
from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes
from orbit.py_linac.lattice_modifications import Replace_Quads_to_OverlappingQuads_Nodes

from orbit.py_linac.lattice import GetGlobalQuadGradient
from orbit.py_linac.lattice import GetGlobalQuadGradientDerivative
from orbit.py_linac.lattice import GetGlobalRF_AxisField

from orbit.py_linac.overlapping_fields import SNS_EngeFunctionFactory

# names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]
names = [
    "DTL6",
]
names = [
    "SCLMed",
]
names = [
    "MEBT",
]

# ---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_linac.xml"

# ---- RF axis fields files location
rf_field_location = "../sns_rf_fields/"

# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.1)


def getQuadFiledsDistribution(seq_name, xml_file_name, rf_field_location, use_overlapped_fields=False, reverse_latt=False):
    accLattice = sns_linac_factory.getLinacAccLattice(
        [
            seq_name,
        ],
        xml_file_name,
    )
    print("debug seq=", seq_name, " Length [m] =", accLattice.getLength())

    if use_overlapped_fields:
        # ---- longitudinal step along the distributed fields lattice
        z_step = 0.005

        # ----------------------------------------------------------------------------------------
        # Here we have three variants:
        # 0. To transformation - hard edge quads and thin RF gaps
        # 1. Replace thin RF-gap node with distributed RF Ez(z) field on the axis
        # 2. Replace hard-edge quads with soft-edge ones, and the quads' fields could overlap
        # 3. Replace RF-gap thin nodes and hard-edge quads with nodes with distributed fields
        # -----------------------------------------------------------------------------------------

        # ---- This function cannot be applied to DTL because the RF fields are overlapping quads
        # if(not seq_name.find("DTL") >= 0):
        # 	Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,rf_field_location,[seq_name,])

        # Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,[seq_name,],[],SNS_EngeFunctionFactory)

        Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(
            accLattice,
            z_step,
            rf_field_location,
            [
                seq_name,
            ],
            [],
            SNS_EngeFunctionFactory,
        )

    # -----------------------------------------------------------------------------
    # Debugging part START
    # -----------------------------------------------------------------------------

    nodes = accLattice.getNodes()
    print("debug =============== DIRECT nodes N=", len(nodes))
    for node_ind in range(len(nodes)):
        node = nodes[node_ind]
        print("debug direct node= %35s " % node.getName(), " pos = %7.3f " % node.getPosition(), " Nchildren=", node.getNumberOfChildren())
        child_nodes = node.getAllChildren()
        for child_node in child_nodes:
            print("                       child = %35s " % child_node.getName())

    # -----------------------------------------------------------------------------
    # Debugging part END
    # -----------------------------------------------------------------------------

    if reverse_latt:
        # ---- reverse the lattice for the backward tracking
        accLattice.reverseOrder()

    # -----------------------------------------------------------------------------
    # Debugging part START
    # -----------------------------------------------------------------------------

    nodes = accLattice.getNodes()
    print("debug =============== REVERSE nodes N=", len(nodes))
    for node_ind in range(len(nodes)):
        node = nodes[node_ind]
        pos = accLattice.getLength() - node.getPosition()
        print("debug reverse node= %35s " % node.getName(), " pos = %7.3f " % pos, " Nchildren=", node.getNumberOfChildren())
        child_nodes = node.getAllChildren()
        for child_node in child_nodes:
            print("                       child = %35s " % child_node.getName())

    # -----------------------------------------------------------------------------
    # Debugging part END
    # -----------------------------------------------------------------------------

    # ---- magn_field_arr[[z,g0,gp0,g1,gp1,Ez],...]
    # ---- g is [T/m] and gp = dq/dz - derivative along z-axis
    z_arr = []
    g_arr = []
    gp_arr = []
    ez_arr = []

    step = 0.001
    n_points = int(accLattice.getLength() / step)
    step = accLattice.getLength() / (n_points - 1)
    for ip in range(n_points):
        z = step * ip
        g = GetGlobalQuadGradient(accLattice, z)
        gp = GetGlobalQuadGradientDerivative(accLattice, z)
        Ez = GetGlobalRF_AxisField(accLattice, z)
        if reverse_latt:
            z = accLattice.getLength() - z
        z_arr.append(z)
        g_arr.append(g)
        gp_arr.append(gp)
        ez_arr.append(Ez)

    return (z_arr, g_arr, gp_arr, ez_arr)


# --------------------------------------------------
#  Plotting the fields
# --------------------------------------------------

import matplotlib.pyplot as plt

# -------------------------------------------------------------------------
# Calculate G[t/m], G'[T/m^2], and Ez(z) for direct and reversed lattices
# -------------------------------------------------------------------------

for seq_name in names:
    (z_arr, g_arr, gp_arr, ez_arr) = getQuadFiledsDistribution(
        seq_name, xml_file_name, rf_field_location, use_overlapped_fields=True, reverse_latt=False
    )

    (z1_arr, g1_arr, gp1_arr, ez1_arr) = getQuadFiledsDistribution(
        seq_name, xml_file_name, rf_field_location, use_overlapped_fields=True, reverse_latt=True
    )

    plt.plot(z_arr, g_arr)
    plt.plot(z1_arr, g1_arr)
    plt.xlabel("z[m]")
    plt.ylabel("G[T/m]")
    # plt.savefig(seq_name+"_g.png")
    plt.show()

    plt.plot(z_arr, gp_arr)
    plt.plot(z1_arr, gp1_arr)
    plt.xlabel("z[m]")
    plt.ylabel("G_prime[T/m^2]")
    plt.savefig(seq_name + "_gp.png")
    plt.show()

    plt.plot(z_arr, ez_arr)
    plt.plot(z1_arr, ez1_arr)
    plt.xlabel("z[m]")
    plt.ylabel("Ez[V/m]")
    # plt.savefig(seq_name+"_ez.png")
    plt.show()


print("Stop.")
