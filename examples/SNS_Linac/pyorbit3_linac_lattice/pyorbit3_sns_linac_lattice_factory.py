#! /usr/bin/env python

"""
This is a test script to check the functionality of the SNS linac parser
and some functions of the linac lattice class
It also will measure how fast we parsing.
"""

import sys
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

# import acc. nodes
from orbit.py_linac.lattice import Quad, Drift, Bend, BaseRF_Gap

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor

# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()

# ---- set the maximal drift length - sometimes it is useful for diagnostics
sns_linac_factory.setMaxDriftLength(0.05)
print("Lattice max drift length =", sns_linac_factory.getMaxDriftLength())

print("=============================================")

# ---- define the sequences in the linac accelerator lattice
names = ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4", "DTL5", "DTL6", "CCL1", "CCL2", "CCL3", "CCL4", "SCLMed", "SCLHigh", "HEBT1", "HEBT2"]

# ---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_linac.xml"
acc_da = XmlDataAdaptor.adaptorForFile(xml_file_name)

# ---- make lattice from XML file - another way to construct a lattice
# acc_lattice_linac = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

acc_lattice_linac = sns_linac_factory.getLinacAccLatticeFromDA(names, acc_da)
seqs = acc_lattice_linac.getSequences()
L_acc = 0.0
for seq in seqs:
    L_acc += seq.getLength()
    print("seq=", seq.getName(), " pos=", seq.getPosition(), " L=", seq.getLength(), " total=", L_acc, " nCav=", len(seq.getRF_Cavities()))

print("=============================================")

cavs = acc_lattice_linac.getRF_Cavities()
print("total number of RF cavities =", len(cavs))

rf_gaps = acc_lattice_linac.getRF_Gaps()
print("total number of RF gaps     =", len(rf_gaps))

quads = acc_lattice_linac.getQuads()
print("total number of Quads       =", len(quads))

print("=============================================")

quads = acc_lattice_linac.getNodesOfClass(Quad)
print("total number of Quads(again)=", len(quads))

quads = acc_lattice_linac.getNodesOfClass(
    Quad,
    [
        "MEBT",
    ],
)
print("total number of Quads (MEBT)=", len(quads))

quads_and_rf_gaps = acc_lattice_linac.getNodesOfClasses(
    [Quad, BaseRF_Gap],
    [
        "MEBT",
    ],
)
print("total number of Quads&Gaps (MEBT) = ", len(quads_and_rf_gaps))


"""
#---- make lattice from XML Data Adaptor and measure the speed
time_start = time.clock()
count = 0
while(1 < 2):
	count += 1
	acc_lattice_linac = sns_linac_factory.getLinacAccLatticeFromDA(names,acc_da)
	tm = (time.clock() - time_start)/count
	print "debug parsing speed =",tm
"""
