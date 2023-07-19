import pytest

"""
This script will track the bunch through MEBT-DTL3 part of
the SNS linac. At the end of MEBT-DTL3 we will add a vacuum window
made of carbon. After tracking through the vacuum window
some partciles will be lost and some deflected.
"""
import os
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

random.seed(100)


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
script_dir = os.path.dirname(os.path.abspath(__file__))
xml_file_name = os.path.join(script_dir, "../sns_linac_xml/sns_linac.xml")
mebt_start_bunch = os.path.join(script_dir, "bunch_at_mebt_start.dat")

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
bunch.readBunch(mebt_start_bunch)

lostbunch = Bunch()
paramsDict = {
    "lostbunch": lostbunch,
}

sns_lattice.trackDesignBunch(bunch)
sns_lattice.trackBunch(bunch, paramsDict)

print("transported bunch size = ", bunch.getSize())
print("       lost bunch size = ", lostbunch.getSize())


def test_bunch():
    assert bunch.getSize() == 864


def test_lostbunch():
    assert lostbunch.getSize() == 136
