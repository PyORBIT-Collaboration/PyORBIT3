#! /usr/bin/env python

"""
This script is an "one-time deal". It was used one time to modify 
the SNS linac XML file. It can be used as example if somebody wants
to do something similar with the XML lattice file.

This is a script will modify the linac lattice file.
It will remove "effLength" parameter from the QUAD type elements, 
because we do not use it in the linac lattice factory.
It also will remove the same parameter from BEND elements and 
specify the "length" parameter equal to "effLength" parameter, 
because we use only "length" in the linac lattice factory.
For DCH and DCV types the "B" parameter was added to "parameters" with
0. value. These are correctors, and B is the default field in [T].
Also the apertures were added to the RF gaps and BENDs.
The aperture parameters are for SNS Linac.
"""

import sys
import time

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor


def setAprtParamsToBEND(node_da):
	if(node_da.hasAttribute("type") and node_da.stringValue("type") == "BEND"):
		params_da = node_da.childAdaptors("parameters")[0]
		params_da.removeAttribute("aperture")
		params_da.setValue("aprt_type",3)
		params_da.setValue("aperture_x",0.420)
		params_da.setValue("aperture_y",0.080)

def setAprtParamsToRF_Gap(accSeq_da,node_da):
	if(node_da.hasAttribute("type") and node_da.stringValue("type") == "RFGAP"):
		params_da = node_da.childAdaptors("parameters")[0]
		if(accSeq_da.getName().find("MEBT") >= 0):
			if(node_da.stringValue("name").find("Bnch01") >= 0 or node_da.stringValue("name").find("Bnch04") >= 0 ):
				params_da.setValue("aprt_type",1)
				params_da.setValue("aperture",0.030)
			if(node_da.stringValue("name").find("Bnch02") >= 0 or node_da.stringValue("name").find("Bnch03") >= 0 ):
				params_da.setValue("aprt_type",1)
				params_da.setValue("aperture",0.036)
		if(accSeq_da.getName().find("DTL") >= 0):	
				params_da.setValue("aprt_type",1)
				params_da.setValue("aperture",0.025)	
		if(accSeq_da.getName().find("CCL") >= 0):	
				params_da.setValue("aprt_type",1)
				params_da.setValue("aperture",0.030)
		if(accSeq_da.getName().find("SCLMed") >= 0):	
				params_da.setValue("aprt_type",1)
				params_da.setValue("aperture",0.086)
		if(accSeq_da.getName().find("SCLHigh") >= 0):	
				params_da.setValue("aprt_type",1)
				params_da.setValue("aperture",0.098)
				

print("==============START=======================")

#---- define the sequences in the linac accelerator lattice
names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]
root_name = "../sns_linac_xml/"

#---- the XML file name with the linac structures
xml_file_names = []
xml_file_names.append("sns_linac.xml")
xml_file_names.append("sns_sts_linac.xml")

for xml_file_name in xml_file_names:
	print("========start file=",xml_file_name)
	acc_da = XmlDataAdaptor.adaptorForFile(root_name+xml_file_name)
	accSeqs_da = acc_da.childAdaptors()
	for accSeq_da in accSeqs_da:
		print("     start seq.=",accSeq_da.getName())
		for node_da in accSeq_da.childAdaptors():
			if(node_da.hasAttribute("type") and (( node_da.stringValue("type") == "QUAD") or ( node_da.stringValue("type") == "BEND"))):
				length = node_da.doubleValue("length")
				params_da = node_da.childAdaptors("parameters")[0]
				if(params_da.hasAttribute("effLength")):
					effLength = params_da.doubleValue("effLength")
					node_da.setValue("length",effLength)
				params_da.removeAttribute("effLength")
				if(xml_file_name.find("_aprt") >= 0): setAprtParamsToBEND(node_da)
			if(node_da.hasAttribute("type") and (( node_da.stringValue("type") == "DCV") or ( node_da.stringValue("type") == "DCH"))):
				params_da = node_da.childAdaptors("parameters")[0]
				params_da.setValue("B",0.)
			if(xml_file_name.find("_aprt") >= 0): setAprtParamsToRF_Gap(accSeq_da,node_da)
	acc_da.writeToFile(root_name+xml_file_name)
		



