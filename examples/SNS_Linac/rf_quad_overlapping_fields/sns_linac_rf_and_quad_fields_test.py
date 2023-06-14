#! /usr/bin/env python

"""
This script will print the data files with quadrupole gradient and
rf fields as functions of the longitudinal positions along each
accelerator sequence in the SNS lattice.
Here we use the SNS specific Enge's functions factory that gives the field
distribution for each quad. User can omit this parameter in the lattice 
transformation calls, and in this case he/she will be using the default factory.
User also can create his/her own factory.

The example includes the three possible lattice transformation functions:
1. Replace_BaseRF_Gap_to_AxisField_Nodes(...) - RF gap replacement
2. Replace_Quads_to_OverlappingQuads_Nodes(...) - Quads replacement
3. Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(...) - RF Gaps and Quads replacement
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

#names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]
names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh"]

#---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_linac.xml"

#---- RF axis fields files location 
dir_location = "../sns_rf_fields/"


#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

#---- output file suffix
output_file_sfx = "_fields_rf_and_quads.dat"

for name in names:
	print("==========================  start of the accSeq=",name)
	#---- make lattice from XML file 
	accLattice = sns_linac_factory.getLinacAccLattice([name,],xml_file_name)
	
	print("Linac initial lattice is ready.  L=",accLattice.getLength())
	
	
	#---- magn_field_arr[[z,g0,gp0,g1,gp1,Ez],...]
	#---- g is [T/m] and gp = dq/dz - derivative along z-axis
	magn_field_arr = []
	step = 0.001
	n_points = int(accLattice.getLength()/step)
	step = accLattice.getLength()/(n_points - 1)
	for ip in range(n_points):
		z = step*ip
		g0 = GetGlobalQuadGradient(accLattice,z)
		gp0 = GetGlobalQuadGradientDerivative(accLattice,z)
		g1 = 0.
		gp1 = 0.
		Ez = 0.
		magn_field_arr.append([z,g0,gp0,g1,gp1,Ez])
	
	#---- longitudinal step along the distributed fields lattice
	z_step = 0.005 
	
	#---- This function cannot be applied to DTL because the RF fields are overlapping quads
	#if(not name.find("DTL") >= 0):
	#	Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,[name,])
	
	#Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,[name,],[],SNS_EngeFunctionFactory)
	
	Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,[name,],[],SNS_EngeFunctionFactory)
		
	print("Linac modified lattice is ready. L=",accLattice.getLength())
	
	for ip in range(n_points):
		[z,g0,gp0,g1,gp1,Ez] = magn_field_arr[ip]
		g1 = GetGlobalQuadGradient(accLattice,z)
		gp1 = GetGlobalQuadGradientDerivative(accLattice,z)
		Ez = GetGlobalRF_AxisField(accLattice,z)
		magn_field_arr[ip] = [z,g0,gp0,g1,gp1,Ez]
	
	fl_out = open(name+output_file_sfx,"w")
	fl_out.write("z[m]        G0[T/m]  dG0/dz[T/m^2]     G1[T/m]   dG1/dz[T/m^2]    Ez[V/m] \n")
	
	for ip in range(n_points):
		[z,g0,gp0,g1,gp1,Ez] = magn_field_arr[ip]
		fl_out.write(" %14.8f   %12.6f  %12.6f     %12.6f  %12.6f    %12.5g "%(z,g0,gp0,g1,gp1,Ez)+"\n")
		
	fl_out.close()
	

print("Stop.")