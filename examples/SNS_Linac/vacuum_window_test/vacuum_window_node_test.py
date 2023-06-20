#! /usr/bin/env python

"""
This script will track the bunch through MEBT-DTL3 part of
the SNS linac. At the end of MEBT-DTL3 we will add a vacuum window
made of carbon. After tracking through the vacuum window
some partciles will be lost and some deflected.
"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor

"""
  material index (0=carbon, 1=aluminum, 2=iron, 3=copper, 4=tantalum, 5=tungstun,
   6=platinum, 7=lead, ma>=8 = black absorber)

  density factor (for materials mixed with air or water). 1.0 for pure

  <accElement length="0.001" name="MEBT:VACWIN" pos="69.5210123" type="VACWIN">
   <parameters material_index="0" density_factor="1.0" />
  </accElement>
"""


def make_vacwin_da(name, length, pos):
    """
    Creates data adaptor for vacuum window node
    """
    node_da = XmlDataAdaptor("accElement")
    node_da.setValue("name", name)
    node_da.setValue("pos", pos)
    node_da.setValue("type", "VACWIN")
    node_da.setValue("length", length)
    params_da = node_da.createChild("parameters")
    params_da.setValue("material_index", 0)
    params_da.setValue("density_factor", 1.0)
    return node_da


def getChildDA(lattice_da, node_name):
    """
    returns the child node of DA with particular name
    """
    for da_ind in range(len(lattice_da.childAdaptors()) - 1):
        node_da = lattice_da.childAdaptors()[da_ind]
        name = node_da.stringValue("name")
        if name == node_name:
            return node_da
    return None


# ==============================================================
#                START of SCRIPT
# ==============================================================
xml_file_name = "../sns_linac_xml/sns_linac.xml"
sns_lattice_da = XmlDataAdaptor.adaptorForFile(xml_file_name)

dtl_da = getChildDA(sns_lattice_da, "DTL3")

# ---- make vacuum window node and add it to the end of the DTL
window_length = 0.005
pos = dtl_da.doubleValue("length") - window_length - 0.00001
vacwin_node = make_vacwin_da("DTL:VACWIN", window_length, pos)
dtl_da.addChildAdaptor(vacwin_node)

# ---- print a structure of the new data adaptor DTL3
# print dtl_da.makeXmlText()

# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.1)

# ---- make lattice from XML file
names = ["MEBT", "DTL1", "DTL2", "DTL3"]
sns_lattice = sns_linac_factory.getLinacAccLatticeFromDA(names, sns_lattice_da)

bunch = Bunch()
bunch.readBunch("bunch_at_mebt_start.dat")

lostbunch = Bunch()
paramsDict = {
    "lostbunch": lostbunch,
}

sns_lattice.trackDesignBunch(bunch)
sns_lattice.trackBunch(bunch, paramsDict)

print("========= Bunch Sizes after Vacuum Window ==========")
print("transported bunch size = ", bunch.getSize())
print("       lost bunch size = ", lostbunch.getSize())
print("====================================================")

bunch.dumpBunch("bunch_at_dtl_exit.dat")
lostbunch.dumpBunch("lostbunch.dat")

print("Stop.")
