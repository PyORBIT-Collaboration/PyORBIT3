#! /usr/bin/env python

"""
This script will track the bunch through the SNS Linac with
the modified lattice. The quads and RF gaps will be replaced by
elements with distributed fields.
"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

# from linac import the C++ RF gap classes
from orbit.core.linac import BaseRfGap, MatrixRfGap, RfGapTTF

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D

from orbit.core.bunch import Bunch, BunchTwissAnalysis

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

# from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice
# from orbit.py_linac.lattice_modifications import Add_rfgap_apertures_to_lattice
# from orbit.py_linac.lattice_modifications import AddMEBTChopperPlatesAperturesToSNS_Lattice
# from orbit.py_linac.lattice_modifications import AddScrapersAperturesToLattice

from orbit.py_linac.lattice import AxisField_and_Quad_RF_Gap
from orbit.py_linac.lattice import OverlappingQuadsNode

from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_to_AxisField_Nodes
from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes
from orbit.py_linac.lattice_modifications import Replace_Quads_to_OverlappingQuads_Nodes

from orbit.py_linac.overlapping_fields import SNS_EngeFunctionFactory

# we take a SNS Linac Bunch generator from a neighboring directory
sys.path.append("../pyorbit3_linac_model")
from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

random.seed(100)

# names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]
# names = ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4", "DTL5", "DTL6", "CCL1", "CCL2", "CCL3", "CCL4", "SCLMed", "SCLHigh", "HEBT1"]
names = ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4"]
# names = ["MEBT",]

# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

# ---- the XML file name with the structure
xml_file_name = "../ess_linac_xml/ess_wdtl_linac.xml"

# ---- make lattice from XML file
accLattice = sns_linac_factory.getLinacAccLattice(names, xml_file_name)

print("Linac lattice is ready. L=", accLattice.getLength())

# ----set up RF Gap Model -------------
# ---- There are three available models at this moment
# ---- BaseRfGap  uses only E0TL*cos(phi)*J0(kr) with E0TL = const
# ---- MatrixRfGap uses a matrix approach like envelope codes
# ---- RfGapTTF uses Transit Time Factors (TTF) like PARMILA
cppGapModel = BaseRfGap
# cppGapModel = MatrixRfGap
# cppGapModel = RfGapTTF
rf_gaps = accLattice.getRF_Gaps()
for rf_gap in rf_gaps:
    rf_gap.setCppGapModel(cppGapModel())

# # ---- If you want to switch off all cavities - remove comments marks
# cavs = accLattice.getRF_Cavities()
# for cav in cavs:
#     cav.setAmp(0)

# ------------------------------------------------------------------
# ---- BaseRF_Gap and Quads will be replaced for specified sequences
# ------------------------------------------------------------------

# ---- longitudinal step along the distributed fields lattice
z_step = 0.005

# ---- axis fields files location
dir_location = "../ess_rf_fields/"

# Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,["MEBT",])

# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6"],[],SNS_EngeFunctionFactory)
# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT","DTL1","DTL2","DTL3",],[],SNS_EngeFunctionFactory)
# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["SCLMed","SCLHigh"],[],SNS_EngeFunctionFactory)
# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT",],[],SNS_EngeFunctionFactory)

# Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,["MEBT",],[],SNS_EngeFunctionFactory)
# Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,["MEBT","DTL1"],[],SNS_EngeFunctionFactory)
# Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,["DTL1",],[],SNS_EngeFunctionFactory)

# ------------Set the linac specific trackers for drifts, quads and RF gaps ------
# ------------It is for situation when the bunch has very big energy spread
# ------------which is unusual
# accLattice.setLinacTracker(True)


# ------------ Add tracking through the longitudinal field component of the quad
# ------------ The longitudinal component is not zero only for the distributed
# ------------ magnetic field of the quad. It means you have to modify lattice with
# ------------ Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes
# ------------ or
# ------------ Replace_Quads_to_OverlappingQuads_Nodes
nodes = accLattice.getNodes()
for node in nodes:
    if isinstance(node, OverlappingQuadsNode) or isinstance(node, AxisField_and_Quad_RF_Gap):
        # node.setUseLongitudinalFieldOfQuad(True)
        pass


# -----------------------------------------------------
# Set up Space Charge Acc Nodes
# -----------------------------------------------------
from orbit.space_charge.sc3d import setSC3DAccNodes, setUniformEllipsesSCAccNodes
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse, SpaceChargeCalc3D

sc_path_length_min = 0.02

print("Set up Space Charge nodes. ")

"""
# set of uniformly charged ellipses Space Charge
nEllipses = 1
calcUnifEllips = SpaceChargeCalcUnifEllipse(nEllipses)
space_charge_nodes = setUniformEllipsesSCAccNodes(accLattice, sc_path_length_min, calcUnifEllips)

"""
# set FFT 3D Space Charge instead of UniformEllipses
sizeX = 16
sizeY = 16
sizeZ = 16
calc3d = SpaceChargeCalc3D(sizeX, sizeY, sizeZ)
space_charge_nodes = setSC3DAccNodes(accLattice, sc_path_length_min, calc3d)

max_sc_length = 0.0
min_sc_length = accLattice.getLength()
for sc_node in space_charge_nodes:
    scL = sc_node.getLengthOfSC()
    # if(scL > 0.028):
    # 	print "debug sc node = ",sc_node.getName()," L=",scL
    if scL > max_sc_length:
        max_sc_length = scL
    if scL < min_sc_length:
        min_sc_length = scL
print("maximal SC length =", max_sc_length, "  min=", min_sc_length)


# print("===== Aperture Nodes START  =======")
# aprtNodes = Add_quad_apertures_to_lattice(accLattice)
# aprtNodes = Add_rfgap_apertures_to_lattice(accLattice,aprtNodes)
# aprtNodes = AddMEBTChopperPlatesAperturesToSNS_Lattice(accLattice,aprtNodes)

# x_size = 0.042
# y_size = 0.042
# aprtNodes = AddScrapersAperturesToLattice(accLattice,"MEBT_Diag:H_SCRP",x_size,y_size,aprtNodes)

# x_size = 0.042
# y_size = 0.042
# aprtNodes = AddScrapersAperturesToLattice(accLattice,"MEBT_Diag:V_SCRP",x_size,y_size,aprtNodes)


# for node in aprtNodes:
# 	print "aprt=",node.getName()," pos =",node.getPosition()

# print("===== Aperture Nodes Added ======= N total=", len(aprtNodes))


# -----TWISS Parameters at the entrance of MEBT ---------------
# transverse emittances are unnormalized and in pi*mm*mrad
# longitudinal emittance is in pi*eV*sec
e_kin_ini = 0.00362  # in [MeV]
mass = 0.938272  # in [GeV]
gamma = (mass + e_kin_ini) / mass
beta = math.sqrt(gamma * gamma - 1.0) / gamma
print("relat. gamma=", gamma)
print("relat.  beta=", beta)
frequency = 352.21e6
v_light = 2.99792458e8  # in [m/sec]

# ------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta
(alphaX, betaX, emittX) = (-0.051805615, 0.20954703, 2.875304356347471)
(alphaY, betaY, emittY) = (-0.30984478, 0.37074849, 2.861819260781874)
(alphaZ, betaZ, emittZ) = (-0.48130325, 0.92564505, 4.079614268808417)

# alphaZ = -alphaZ

# ---make emittances un-normalized XAL units [m*rad]
emittX = 1.0e-6 * emittX  # / (gamma * beta)
emittY = 1.0e-6 * emittY  # / (gamma * beta)
emittZ = 1.0e-6 * emittZ  # / (gamma**3 * beta)
print(" ========= XAL Twiss ===========")
print(" alpha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f " % (alphaX, betaX, emittX * 1.0e6))
print(" alpha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f " % (alphaY, betaY, emittY * 1.0e6))
print(" alpha beta emitt[mm*mrad] Z= %6.4f %6.4f %6.4f " % (alphaZ, betaZ, emittZ * 1.0e6))

# ---- long. size in mm
sizeZ = math.sqrt(emittZ * betaZ) * 1.0e3

# ---- transform to pyORBIT emittance[GeV*m]
emittZ = emittZ * gamma**3 * beta**2 * mass
betaZ = betaZ / (gamma**3 * beta**2 * mass)

print(" ========= PyORBIT Twiss ===========")
print(" alpha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f " % (alphaX, betaX, emittX * 1.0e6))
print(" alpha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f " % (alphaY, betaY, emittY * 1.0e6))
print(" alpha beta emitt[mm*MeV] Z= %6.4f %6.4f %6.4f " % (alphaZ, betaZ, emittZ * 1.0e6))

twissX = TwissContainer(alphaX, betaX, emittX)
twissY = TwissContainer(alphaY, betaY, emittY)
twissZ = TwissContainer(alphaZ, betaZ, emittZ)

print("Start Bunch Generation.")
bunch_gen = SNS_Linac_BunchGenerator(twissX, twissY, twissZ, frequency=frequency)
bunch_gen.bunch.mass(mass)
bunch_gen.bunch.charge(1.0)

# set the initial kinetic energy in GeV
bunch_gen.setKinEnergy(e_kin_ini)

# set the beam peak current in mA
bunch_gen.setBeamCurrent(62.5)

# bunch_in = bunch_gen.getBunch(nParticles=100000, distributorClass=WaterBagDist3D)
bunch_in = bunch_gen.getBunch(nParticles=100000, distributorClass=GaussDist3D)
# bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = KVDist3D)

print("Bunch Generation completed.")

# set up design
accLattice.trackDesignBunch(bunch_in)

print("Design tracking completed.")

# track through the lattice
paramsDict = {"old_pos": -1.0, "count": 0, "pos_step": 0.01}
actionContainer = AccActionsContainer("Bunch Tracking")

pos_start = 0.0

twiss_analysis = BunchTwissAnalysis()

file_out = open("pyorbit_twiss_sizes_and_ekin.dat", "w")

s = " Node   position "
s += "   alphaX betaX emittX  normEmittX"
s += "   alphaY betaY emittY  normEmittY"
s += "   alphaZ betaZ emittZ  emittZphiMeV"
s += "   sizeX sizeY sizeZ_deg"
s += "   eKin Nparts "
file_out.write(s + "\n")
print(" N node   position    sizeX  sizeY  sizeZdeg  eKin Nparts ")


def action_entrance(paramsDict):
    node = paramsDict["node"]
    bunch = paramsDict["bunch"]
    pos = paramsDict["path_length"]
    if paramsDict["old_pos"] == pos:
        return
    if paramsDict["old_pos"] + paramsDict["pos_step"] > pos:
        return
    paramsDict["old_pos"] = pos
    paramsDict["count"] += 1
    gamma = bunch.getSyncParticle().gamma()
    beta = bunch.getSyncParticle().beta()
    twiss_analysis.analyzeBunch(bunch)
    x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1] * twiss_analysis.getTwiss(0)[3]) * 1000.0
    y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1] * twiss_analysis.getTwiss(1)[3]) * 1000.0
    z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1] * twiss_analysis.getTwiss(2)[3]) * 1000.0
    z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
    z_rms_deg = z_to_phase_coeff * z_rms / 1000.0
    nParts = bunch.getSizeGlobal()
    (alphaX, betaX, emittX) = (twiss_analysis.getTwiss(0)[0], twiss_analysis.getTwiss(0)[1], twiss_analysis.getTwiss(0)[3] * 1.0e6)
    (alphaY, betaY, emittY) = (twiss_analysis.getTwiss(1)[0], twiss_analysis.getTwiss(1)[1], twiss_analysis.getTwiss(1)[3] * 1.0e6)
    (alphaZ, betaZ, emittZ) = (twiss_analysis.getTwiss(2)[0], twiss_analysis.getTwiss(2)[1], twiss_analysis.getTwiss(2)[3] * 1.0e6)
    norm_emittX = emittX * gamma * beta
    norm_emittY = emittY * gamma * beta
    # ---- phi_de_emittZ will be in [pi*deg*MeV]
    phi_de_emittZ = z_to_phase_coeff * emittZ
    eKin = bunch.getSyncParticle().kinEnergy() * 1.0e3
    s = " %35s  %4.5f " % (node.getName(), pos + pos_start)
    s += "   %6.4f  %6.4f  %6.4f  %6.4f   " % (alphaX, betaX, emittX, norm_emittX)
    s += "   %6.4f  %6.4f  %6.4f  %6.4f   " % (alphaY, betaY, emittY, norm_emittY)
    s += "   %6.4f  %6.4f  %6.4f  %6.4f   " % (alphaZ, betaZ, emittZ, phi_de_emittZ)
    s += "   %5.3f  %5.3f  %5.3f " % (x_rms, y_rms, z_rms_deg)
    s += "  %10.6f   %8d " % (eKin, nParts)
    file_out.write(s + "\n")
    file_out.flush()
    s_prt = " %5d  %35s  %4.5f " % (paramsDict["count"], node.getName(), pos + pos_start)
    s_prt += "  %5.3f  %5.3f   %5.3f " % (x_rms, y_rms, z_rms_deg)
    s_prt += "  %10.6f   %8d " % (eKin, nParts)
    print(s_prt)


def action_exit(paramsDict):
    action_entrance(paramsDict)


actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

time_start = time.process_time()

results_arr = accLattice.trackBunch(bunch_in, paramsDict=paramsDict, actionContainer=actionContainer)

time_exec = time.process_time() - time_start
print("time[sec]=", time_exec)

# ---- dump the final bunch
# bunch_in.dumpBunch("bunch_test.dat")

file_out.close()
