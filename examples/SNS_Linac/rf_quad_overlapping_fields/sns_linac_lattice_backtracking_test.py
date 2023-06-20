#! /usr/bin/env python

"""
This script will track back and forth the bunch through the SNS Linac with
the modified lattice. The quads and RF gaps will be replaced by
elements with distributed fields. The results of tracking should be
reversible.
"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

# from linac import the C++ RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D

from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice
from orbit.py_linac.lattice_modifications import Add_rfgap_apertures_to_lattice
from orbit.py_linac.lattice_modifications import AddMEBTChopperPlatesAperturesToSNS_Lattice
from orbit.py_linac.lattice_modifications import AddScrapersAperturesToLattice

from orbit.py_linac.lattice import AxisField_and_Quad_RF_Gap
from orbit.py_linac.lattice import OverlappingQuadsNode

from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_to_AxisField_Nodes
from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes
from orbit.py_linac.lattice_modifications import Replace_Quads_to_OverlappingQuads_Nodes

from orbit.py_linac.overlapping_fields import SNS_EngeFunctionFactory

# we take a SNS Linac Bunch generator from a neighboring directory
sys.path.append("../pyorbit_linac_model")
from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

from bunch import Bunch

random.seed(100)

# names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]
# names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1"]
# names = ["MEBT","DTL1","DTL2","DTL3",]
# names = ["MEBT","DTL1"]
# names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed"]
names = ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4", "DTL5", "DTL6", "CCL1", "CCL2", "CCL3", "CCL4", "SCLMed", "SCLHigh"]
# names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6"]
# names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4"]

# names = ["MEBT",]

# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

# ---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_linac.xml"

# ---- make lattice from XML file
accLattice = sns_linac_factory.getLinacAccLattice(names, xml_file_name)

print("Linac lattice is ready. L=", accLattice.getLength())

print("===========Cavities' and gaps' pphases =================")

cavs = accLattice.getRF_Cavities()
for cav in cavs:
    cav_phase = cav.getPhase()
    rf_gaps = cav.getRF_GapNodes()
    print("debug cav = ", cav.getName(), " phase [deg] =", cav_phase * 180.0 / math.pi)
    for rf_gap in rf_gaps:
        print("      debug rf_gap = ", rf_gap.getName(), " gap phase [deg] =", rf_gap.getGapPhase() * 180.0 / math.pi)

print("=============================================")

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

# ---- If you want to switch off all cavities - remove comments marks
# cavs = accLattice.getRF_Cavities()
# for cav in cavs:
# 	cav.setAmp(0.)

# ------------------------------------------------------------------
# ---- BaseRF_Gap and Quads will be replaced for specified sequences
# ------------------------------------------------------------------

# ---- longitudinal step along the distributed fields lattice
z_step = 0.005

# ---- axis fields files location
dir_location = "../sns_rf_fields/"

# Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,["MEBT",])

# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6"],[],SNS_EngeFunctionFactory)
# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT","DTL1","DTL2","DTL3",],[],SNS_EngeFunctionFactory)
# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["SCLMed","SCLHigh"],[],SNS_EngeFunctionFactory)
# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT",],[],SNS_EngeFunctionFactory)
# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT","DTL1"],[],SNS_EngeFunctionFactory)
# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,names,[],SNS_EngeFunctionFactory)


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
from spacecharge import SpaceChargeCalcUnifEllipse, SpaceChargeCalc3D

sc_path_length_min = 0.02

print("Set up Space Charge nodes. ")

# set of uniformly charged ellipses Space Charge
nEllipses = 1
calcUnifEllips = SpaceChargeCalcUnifEllipse(nEllipses)
space_charge_nodes = setUniformEllipsesSCAccNodes(accLattice, sc_path_length_min, calcUnifEllips)

"""
# set FFT 3D Space Charge instead of UniformEllipses
sizeX = 64
sizeY = 64
sizeZ = 64
calc3d = SpaceChargeCalc3D(sizeX,sizeY,sizeZ)
space_charge_nodes =  setSC3DAccNodes(accLattice,sc_path_length_min,calc3d)
"""

max_sc_length = 0.0
min_sc_length = accLattice.getLength()
for sc_node in space_charge_nodes:
    scL = sc_node.getLengthOfSC()
    if scL > max_sc_length:
        max_sc_length = scL
    if scL < min_sc_length:
        min_sc_length = scL
print("maximal SC length =", max_sc_length, "  min=", min_sc_length)


print("===== Aperture Nodes START  =======")
aprtNodes = Add_quad_apertures_to_lattice(accLattice)
# aprtNodes = Add_rfgap_apertures_to_lattice(accLattice,aprtNodes)
# aprtNodes = AddMEBTChopperPlatesAperturesToSNS_Lattice(accLattice,aprtNodes)

x_size = 0.042
y_size = 0.042
# aprtNodes = AddScrapersAperturesToLattice(accLattice,"MEBT_Diag:H_SCRP",x_size,y_size,aprtNodes)

x_size = 0.042
y_size = 0.042
# aprtNodes = AddScrapersAperturesToLattice(accLattice,"MEBT_Diag:V_SCRP",x_size,y_size,aprtNodes)


# for node in aprtNodes:
# 	print "aprt=",node.getName()," pos =",node.getPosition()

print("===== Aperture Nodes Added ======= N total=", len(aprtNodes))


# -----TWISS Parameters at the entrance of MEBT ---------------
# transverse emittances are unnormalized and in pi*mm*mrad
# longitudinal emittance is in pi*eV*sec
e_kin_ini = 0.0025  # in [GeV]
mass = 0.939294  # in [GeV]
gamma = (mass + e_kin_ini) / mass
beta = math.sqrt(gamma * gamma - 1.0) / gamma
print("relat. gamma=", gamma)
print("relat.  beta=", beta)
frequency = 402.5e6
v_light = 2.99792458e8  # in [m/sec]

# ------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta
(alphaX, betaX, emittX) = (-1.9620, 0.1831, 0.21)
(alphaY, betaY, emittY) = (1.7681, 0.1620, 0.21)
(alphaZ, betaZ, emittZ) = (0.0196, 0.5844, 0.24153)

alphaZ = -alphaZ

# ---make emittances un-normalized XAL units [m*rad]
emittX = 1.0e-6 * emittX / (gamma * beta)
emittY = 1.0e-6 * emittY / (gamma * beta)
emittZ = 1.0e-6 * emittZ / (gamma**3 * beta)
print(" ========= XAL Twiss ===========")
print(" aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f " % (alphaX, betaX, emittX * 1.0e6))
print(" aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f " % (alphaY, betaY, emittY * 1.0e6))
print(" aplha beta emitt[mm*mrad] Z= %6.4f %6.4f %6.4f " % (alphaZ, betaZ, emittZ * 1.0e6))

# ---- long. size in mm
sizeZ = math.sqrt(emittZ * betaZ) * 1.0e3

# ---- transform to pyORBIT emittance[GeV*m]
emittZ = emittZ * gamma**3 * beta**2 * mass
betaZ = betaZ / (gamma**3 * beta**2 * mass)

print(" ========= PyORBIT Twiss ===========")
print(" aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f " % (alphaX, betaX, emittX * 1.0e6))
print(" aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f " % (alphaY, betaY, emittY * 1.0e6))
print(" aplha beta emitt[mm*MeV] Z= %6.4f %6.4f %6.4f " % (alphaZ, betaZ, emittZ * 1.0e6))

twissX = TwissContainer(alphaX, betaX, emittX)
twissY = TwissContainer(alphaY, betaY, emittY)
twissZ = TwissContainer(alphaZ, betaZ, emittZ)

print("Start Bunch Generation.")
bunch_gen = SNS_Linac_BunchGenerator(twissX, twissY, twissZ)

# set the initial kinetic energy in GeV
bunch_gen.setKinEnergy(e_kin_ini)

# set the beam peak current in mA
bunch_gen.setBeamCurrent(38.0)

bunch_rfq = bunch_gen.getBunch(nParticles=5000, distributorClass=WaterBagDist3D)
# bunch_rfq = bunch_gen.getBunch(nParticles = 5000, distributorClass = GaussDist3D)
# bunch_rfq = bunch_gen.getBunch(nParticles = 5000, distributorClass = KVDist3D)

print("Bunch Generation completed.")

bunch_forward = Bunch()
bunch_rfq.copyBunchTo(bunch_forward)


# ----------------------------------------------
#  Bunch tracking function
# ----------------------------------------------
def TrackingBunch(accLattice, bunch, print_info=False):
    # set up design
    accLattice.trackDesignBunch(bunch)

    # track through the lattice
    paramsDict = {"old_pos": -1.0, "count": 0, "pos_step": 0.01}
    actionContainer = AccActionsContainer("Bunch Tracking")

    pos_start = 0.0

    twiss_analysis = BunchTwissAnalysis()

    results_arr = []

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
        s_prt = " %5d  %35s  %4.5f " % (paramsDict["count"], node.getName(), pos + pos_start)
        s_prt += "  %5.3f  %5.3f   %5.3f " % (x_rms, y_rms, z_rms_deg)
        s_prt += "  %10.6f   %8d " % (eKin, nParts)
        if print_info:
            print(s_prt)
        twiss_arr = [(alphaX, betaX, emittX, norm_emittX), (alphaY, betaY, emittY, norm_emittY), (alphaZ, betaZ, emittZ, phi_de_emittZ)]
        rms_arr = [x_rms, y_rms, z_rms_deg]
        results_arr.append([node.getName(), pos, rms_arr, twiss_arr, eKin, nParts])

    def action_exit(paramsDict):
        action_entrance(paramsDict)

    actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
    actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

    accLattice.trackBunch(bunch, paramsDict=paramsDict, actionContainer=actionContainer)

    return results_arr


print("====== Start forward tracking ")

results_arr = TrackingBunch(accLattice, bunch_forward, False)


def BunchTransformerFunc(bunch):
    """
    This function will reverse all xp, yp, z coordinates of the bunch.
    We have to change the sign of the z because the tail will be the head
    of the bunch, but the sign of dE will not change because of the same reason.
    """
    nParts = bunch.getSize()
    for i in range(nParts):
        (xp, yp, z, dE) = (bunch.xp(i), bunch.yp(i), bunch.z(i), bunch.dE(i))
        bunch.xp(i, -xp)
        bunch.yp(i, -yp)
        bunch.z(i, -z)
        # --- dE should not change the sign
        # bunch.dE(i,-dE)
    bunch.getSyncParticle().time(0.0)


BunchTransformerFunc(bunch_forward)
bunch_backward = bunch_forward


# ---- reverse the lattice for the backward tracking
accLattice.reverseOrder()

print("========== The lattice was reversed =================")

cavs = accLattice.getRF_Cavities()
for cav in cavs:
    cav_phase = cav.getPhase()
    rf_gaps = cav.getRF_GapNodes()
    print("debug cav = ", cav.getName(), " phase [deg] =", cav_phase * 180.0 / math.pi)
    for rf_gap in rf_gaps:
        print("      debug rf_gap = ", rf_gap.getName(), " gap phase [deg] =", rf_gap.getGapPhase() * 180.0 / math.pi)

print("=============================================")
print("====== Start backward tracking ")

results_backward_arr = TrackingBunch(accLattice, bunch_backward, False)


def GetPlottingArrays(results_arr, length, order=+1):
    pos_arr = []
    rms_x_arr = []
    rms_y_arr = []
    rms_phi_arr = []
    for [name, pos, rms_arr, twiss_arr, eKin, nParts] in results_arr:
        if order < 0:
            pos = length - pos
        pos_arr.append(pos)
        rms_x_arr.append(rms_arr[0])
        rms_y_arr.append(rms_arr[1])
        rms_phi_arr.append(rms_arr[2])
    return (pos_arr, rms_x_arr, rms_y_arr, rms_phi_arr)


# ---------------------------
# Plot part
# ---------------------------
import matplotlib.pyplot as plt

length = accLattice.getLength()

(pos_arr, rms_x_arr, rms_y_arr, rms_phi_arr) = GetPlottingArrays(results_arr, length, order=+1)
(pos_r_arr, rms_x_r_arr, rms_y_r_arr, rms_phi_r_arr) = GetPlottingArrays(results_backward_arr, length, order=-1)

(line_1,) = plt.plot(pos_arr, rms_x_arr, label="Forward Tracking")
(line_2,) = plt.plot(pos_r_arr, rms_x_r_arr, label="Backward Tracking")
plt.legend(handles=[line_1, line_2])
plt.title("Horizontal RMS Sizes")
plt.ylabel("RMS x, mm")
plt.xlabel("pos,m ")

plt.show()

(line_1,) = plt.plot(pos_arr, rms_y_arr, label="Forward Tracking")
(line_2,) = plt.plot(pos_r_arr, rms_y_r_arr, label="Backward Tracking")
plt.legend(handles=[line_1, line_2])
plt.title("Vertical RMS Sizes")
plt.ylabel("RMS y, mm")
plt.xlabel("pos,m ")

plt.show()

(line_1,) = plt.plot(pos_arr, rms_phi_arr, label="Forward Tracking")
(line_2,) = plt.plot(pos_r_arr, rms_phi_r_arr, label="Backward Tracking")
plt.legend(handles=[line_1, line_2])
plt.title("Longitudinal RMS Sizes")
plt.ylabel("RMS Phi, deg")
plt.xlabel("pos,m ")


plt.show()


print("Stop")
sys.exit(0)
