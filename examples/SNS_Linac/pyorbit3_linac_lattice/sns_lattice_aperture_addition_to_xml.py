#! /usr/bin/env python

"""
This script is an "one-time deal". It was used one time to modify
the SNS linac XML file. It can be used as example if somebody wants
to do something similar with the XML lattice file.

This is a script to add aperture parameters to the SNS linac xml file.
"""

import sys
import time

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor


print("==============START=======================")

# ---- define the sequences in the linac accelerator lattice
names = ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4", "DTL5", "DTL6", "CCL1", "CCL2", "CCL3", "CCL4", "SCLMed", "SCLHigh", "HEBT1", "HEBT2"]

# ---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_linac.xml"
acc_da = XmlDataAdaptor.adaptorForFile(xml_file_name)

# -----------MEBT

quad_apert_dict = {}
quad_apert_dict["MEBT_Mag:QH01"] = "0.032"
quad_apert_dict["MEBT_Mag:QV02"] = "0.032"
quad_apert_dict["MEBT_Mag:QH03"] = "0.032"
quad_apert_dict["MEBT_Mag:QV04"] = "0.032"
quad_apert_dict["MEBT_Mag:QH05"] = "0.042"
quad_apert_dict["MEBT_Mag:QV06"] = "0.042"
quad_apert_dict["MEBT_Mag:QH07"] = "0.042"
quad_apert_dict["MEBT_Mag:QH08"] = "0.042"
quad_apert_dict["MEBT_Mag:QV09"] = "0.042"
quad_apert_dict["MEBT_Mag:QH10"] = "0.042"
quad_apert_dict["MEBT_Mag:QV11"] = "0.032"
quad_apert_dict["MEBT_Mag:QH12"] = "0.032"
quad_apert_dict["MEBT_Mag:QV13"] = "0.032"
quad_apert_dict["MEBT_Mag:QH14"] = "0.032"

names = [
    "MEBT",
]

accSeqNodes = acc_da.childAdaptors()
for accSeqNode in accSeqNodes:
    seqName = accSeqNode.getName()
    if seqName not in names:
        continue
    print("seq=", seqName)
    for node in accSeqNode.childAdaptors():
        if node.hasAttribute("type") and (node.stringValue("type") == "QUAD"):
            params_nodes = node.childAdaptors("parameters")
            if len(params_nodes) == 1:
                param_da = params_nodes[0]
                quad_name = node.stringValue("name")
                aperture_val_str = quad_apert_dict[quad_name]
                param_da.setValue("aperture", aperture_val_str)
                param_da.setValue("aprt_type", "1")
                print("   quad=", quad_name, " d=", aperture_val_str)


# -----------DTL CCL
names = ["DTL1", "DTL2", "DTL3", "DTL4", "DTL5", "DTL6", "CCL1", "CCL2", "CCL3", "CCL4"]
quad_apert_dict = {}
quad_apert_dict["DTL1"] = "0.025"
quad_apert_dict["DTL2"] = "0.025"
quad_apert_dict["DTL3"] = "0.025"
quad_apert_dict["DTL4"] = "0.025"
quad_apert_dict["DTL5"] = "0.025"
quad_apert_dict["DTL6"] = "0.025"

quad_apert_dict["CCL1"] = "0.030"
quad_apert_dict["CCL2"] = "0.030"
quad_apert_dict["CCL3"] = "0.030"
quad_apert_dict["CCL4"] = "0.030"

accSeqNodes = acc_da.childAdaptors()
for accSeqNode in accSeqNodes:
    seqName = accSeqNode.getName()
    if seqName not in names:
        continue
    print("seq=", seqName)
    for node in accSeqNode.childAdaptors():
        if node.hasAttribute("type") and (node.stringValue("type") == "QUAD"):
            params_nodes = node.childAdaptors("parameters")
            if len(params_nodes) == 1:
                param_da = params_nodes[0]
                quad_name = node.stringValue("name")
                aperture_val_str = quad_apert_dict[seqName]
                param_da.setValue("aperture", aperture_val_str)
                param_da.setValue("aprt_type", "1")
                print("   quad=", quad_name, " d=", aperture_val_str)

# ----------- SCL
names = ["SCLMed", "SCLHigh"]

accSeqNodes = acc_da.childAdaptors()
for accSeqNode in accSeqNodes:
    seqName = accSeqNode.getName()
    if seqName not in names:
        continue
    print("seq=", seqName)
    for node in accSeqNode.childAdaptors():
        if node.hasAttribute("type") and (node.stringValue("type") == "QUAD"):
            params_nodes = node.childAdaptors("parameters")
            if len(params_nodes) == 1:
                param_da = params_nodes[0]
                quad_name = node.stringValue("name")
                aperture_val_str = "0.080"
                if quad_name == "SCL_Mag:QH00" or quad_name == "SCL_Mag:QV00":
                    aperture_val_str = "0.048"
                param_da.setValue("aperture", aperture_val_str)
                param_da.setValue("aprt_type", "1")
                print("   quad=", quad_name, " d=", aperture_val_str)

# -----------HEBT1 HEBT2
names = ["HEBT1", "HEBT2"]
quad_apert_dict = {}
quad_apert_dict["HEBT1"] = "0.240"
quad_apert_dict["HEBT2"] = "0.420"

accSeqNodes = acc_da.childAdaptors()
for accSeqNode in accSeqNodes:
    seqName = accSeqNode.getName()
    if seqName not in names:
        continue
    print("seq=", seqName)
    for node in accSeqNode.childAdaptors():
        if node.hasAttribute("type") and (node.stringValue("type") == "QUAD"):
            params_nodes = node.childAdaptors("parameters")
            if len(params_nodes) == 1:
                param_da = params_nodes[0]
                quad_name = node.stringValue("name")
                aperture_val_str = quad_apert_dict[seqName]
                param_da.setValue("aperture", aperture_val_str)
                param_da.setValue("aprt_type", "1")
                print("   quad=", quad_name, " d=", aperture_val_str)

# ---------DUMP XML file
acc_da.writeToFile("../sns_linac_xml/sns_linac_with_aprt.xml")

print("==============STOP=======================")
