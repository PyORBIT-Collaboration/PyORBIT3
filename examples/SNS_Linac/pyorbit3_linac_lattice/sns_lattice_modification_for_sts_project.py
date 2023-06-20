#! /usr/bin/env python

"""
This script is an "one-time deal". It was used one time to modify
the SNS linac XML file. It can be used as example if somebody wants
to do something similar with the XML lattice file.

This is a script to modify the linac lattice file for the Second Target Station
project. It will add additional SCL cavities for modules 23-32.
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

accSeq_da = acc_da.childAdaptors("SCLHigh")[0]

dcv_da_arr = []
count = 11
pos_old = 0.0
for node in accSeq_da.childAdaptors():
    if node.hasAttribute("type") and (node.stringValue("type") == "DCV"):
        pos = node.doubleValue("pos")
        delta_pos = pos - pos_old
        pos_old = pos
        count += 1
        if count > 23:
            dcv_da_arr.append(node)
            print("node=", node.stringValue("name"), " pos = %12.5f " % pos, " delta=  %12.5f " % delta_pos)


def getIndexAndNode(da, name):
    arr = da.childAdaptors()
    for node in arr:
        if node.stringValue("name") == name:
            index = arr.index(node)
            return (index, node)


(ind_dcv22, dcv22_node_da) = getIndexAndNode(accSeq_da, "SCL_Mag:DCV22")
pos_0 = accSeq_da.childAdaptors()[ind_dcv22].doubleValue("pos")
print("ind=", ind_dcv22, " node=", accSeq_da.childAdaptors()[ind_dcv22].stringValue("name"))

gaps_arr = []
for ind in range(ind_dcv22 + 1, ind_dcv22 + 24 + 1):
    node = accSeq_da.childAdaptors()[ind]
    pos = node.doubleValue("pos")
    pos_new = pos - pos_0
    print("ind=", ind, " node=", node.stringValue("name"), " pos = %12.5f " % node.doubleValue("pos"), " delta= %12.5f" % pos_new)
    node_new = node.getDeepCopy()
    node_new.setValue("pos", pos_new)
    gaps_arr.append(node_new)

cav_new_names = []
for ind in range(23, 32):
    if ind == 23 or ind == 25:
        continue
    ind_mod = ind + 1
    (ind_dcv, dcv_node_da) = getIndexAndNode(accSeq_da, "SCL_Mag:DCV" + str(ind))
    # print "debug ind=",ind," dcv name=",dcv_node_da.stringValue("name")," pos = %12.5f "%dcv_node_da.doubleValue("pos")," ind_mod=",ind_mod," cav=","SCL:Cav"+str(ind_mod)+"a"
    pos0 = dcv_node_da.doubleValue("pos")
    letter_arr = ["a", "b", "c", "d"]
    gap_count = 0
    for ind_cav in range(4):
        letter = letter_arr[ind_cav]
        cav_name = "SCL:Cav" + str(ind_mod) + letter
        print("debug =================================== new cavity =", cav_name)
        pos_cav_avg = 0.0
        for ind_rf_gap in range(6):
            gap_da = gaps_arr[gap_count].getDeepCopy()
            gap_name = "SCL_RF:Cav" + str(ind_mod) + letter + ":" + gap_da.stringValue("name").split(":")[2]
            gap_da.setValue("name", gap_name)
            pos = pos0 + gap_da.doubleValue("pos")
            pos_cav_avg += pos
            gap_da.setValue("pos", float("%12.5f" % pos))
            gap_da.childAdaptors("parameters")[0].setValue("cavity", cav_name)
            accSeq_da.childAdaptors().insert(gap_count + ind_dcv + 1, gap_da)
            print("debug cav=", cav_name, " rf gap=", gap_da.stringValue("name"), " pos=", gap_da.doubleValue("pos"))
            gap_count += 1
        cav_new_names.append([cav_name, float("%12.5f" % (pos_cav_avg / 6.0))])

count_mod = 0
for ind in range(len(cav_new_names)):
    if ind % 4 == 0:
        count_mod += 1
        print("====== new mod. index =", count_mod)
    [cav_name, pos_avg] = cav_new_names[ind]
    print("cav=", cav_name, " pos=", pos_avg)

# ---- we will set all amplitudes for cavities to 0 for now
cavs_da = accSeq_da.childAdaptors("Cavities")[0]
for ind in range(len(cav_new_names)):
    [cav_name, pos_avg] = cav_new_names[ind]
    cav_da = cavs_da.childAdaptors()[0].getDeepCopy()
    cav_da.setValue("name", cav_name)
    cav_da.setValue("pos", pos_avg)
    cav_da.setValue("ampl", 0.0)
    cavs_da.addChildAdaptor(cav_da)

pos_old = 0.0
for cav_da in cavs_da.childAdaptors():
    pos = cav_da.doubleValue("pos")
    delta = pos - pos_old
    print("debug cav=", cav_da.stringValue("name"), " delta=", delta)
    pos_old = pos


acc_da.writeToFile("../sns_linac_xml/sns_sts_linac.xml")


print("Stop.")
