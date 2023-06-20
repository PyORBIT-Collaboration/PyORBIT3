#! /usr/bin/env python

"""
This script is an "one-time deal". It was used one time to modify
the SNS linac and STS XML file. It can be used as example if somebody wants
to do something similar with the XML lattice file.

This is a script will add the MARKER lattice elements at the beginning and at the end
of each sequence in the lattice.
The name of the elements will be <sequece name>_START and <sequece name>_END
"""

import sys
import time

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor


def getFirstInd(nodes, name):
    ind = 0
    for node in nodes:
        if node.getName() == name:
            return ind
        ind += 1
    return -1


print("==============START SCRIPT=======================")

# ---- define the sequences in the linac accelerator lattice
names = ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4", "DTL5", "DTL6", "CCL1", "CCL2", "CCL3", "CCL4", "SCLMed", "SCLHigh", "HEBT1", "HEBT2"]

# ---- the XML file name with the structure
xml_file_names = ["../sns_linac_xml/sns_linac.xml", "../sns_linac_xml/sns_sts_linac.xml"]
for xml_file_name in xml_file_names:
    print("file =", xml_file_name)
    acc_da = XmlDataAdaptor.adaptorForFile(xml_file_name)

    for name in names:
        accSeq_da = acc_da.childAdaptors(name)[0]
        nodes = accSeq_da.childAdaptors()
        # ==================
        start_marker = XmlDataAdaptor("accElement")
        start_marker.setValue("name", name + "_START")
        start_marker.setValue("length", "0.0")
        start_marker.setValue("pos", "0.0")
        start_marker.setValue("type", "MARKER")
        param_da = XmlDataAdaptor("parameters")
        start_marker.addChildAdaptor(param_da)
        # ===================
        end_marker = XmlDataAdaptor("accElement")
        end_marker.setValue("name", name + "_END")
        end_marker.setValue("length", "0.0")
        end_marker.setValue("pos", accSeq_da.stringValue("length"))
        end_marker.setValue("type", "MARKER")
        param_da = XmlDataAdaptor("parameters")
        end_marker.addChildAdaptor(param_da)
        # ===================
        ind_end = getFirstInd(nodes, "Cavities")
        if ind_end >= 0:
            accSeq_da.addChildAdaptor(end_marker, ind_end)
        else:
            accSeq_da.addChildAdaptor(end_marker)
        accSeq_da.addChildAdaptor(start_marker, 0)

    acc_da.writeToFile(xml_file_name)

print("Stop.")
