#! /usr/bin/env python

"""
This script will track the bunch through the SNS Linac with an upgraded
for the second target station (STS) SCL linac.

SCL Second Target Station linac does not have cavities in the SCL modules
24 and 26.

The model includes the following features:
1. At the beginning you have to generate the bunch at SCL:23d:Rg06 exit by
   un-commenting the part with 'track bunch to the last RF gap in SCL:Cav23d'.
   The script will generate the bunch and will stop. Then comment this part back.
   The rest of the script will read the bunch from the file and will track it
   to end of HEBT2.

2. User will control which cavities are included in acceleration by changing the
   parts of the script with comments 'definitions of the activated cavities with Amp = 1'
   By default the cavities after SCL:23d have Amp=0, and they are not changing the beam.

3. At the end the script will plot the longitudinal phase space by using
   the Gnuplot package.

4. The script includes apertures for all quads

5. There are no overlapping fields quads in the MEBT.

"""

import sys
import math
import random
import time
import orbit.core

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

# from linac import the C++ RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D

from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice
from orbit.py_linac.lattice_modifications import Add_rfgap_apertures_to_lattice
from orbit.py_linac.lattice_modifications import AddMEBTChopperPlatesAperturesToSNS_Lattice
from orbit.py_linac.lattice_modifications import AddScrapersAperturesToLattice

from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from sns_linac_bunch_generator import SNS_Linac_BunchGenerator

# import the utilities
from orbit.utils import phaseNearTargetPhase, phaseNearTargetPhaseDeg


def setSynchPhase(bunch_in, accLattice, cav_name, synchPhaseDeg):
    """
    This function will find the first RF gap phase to get the average
    phase for all gaps equal to the specified synchronous phase.
    """
    # print "debug start def setSynchPhase(...)"
    b = Bunch()
    rf_cav = accLattice.getRF_Cavity(cav_name)
    rf_gaps = rf_cav.getRF_GapNodes()
    ind_start = accLattice.getNodeIndex(rf_gaps[0])
    ind_stop = accLattice.getNodeIndex(rf_gaps[len(rf_gaps) - 1])
    bunch_in.copyEmptyBunchTo(b)
    b = accLattice.trackDesignBunch(b, None, None, -1, ind_start - 1)
    e_kin_in = b.getSyncParticle().kinEnergy()
    e_kin_max = 0
    phase_max = 0.0
    for ind in range(-180, 180):
        phase_deg = ind * 1.0
        rf_cav.setPhase(phase_deg * math.pi / 180.0)
        bunch_in.copyEmptyBunchTo(b)
        b.getSyncParticle().kinEnergy(e_kin_in)
        b = accLattice.trackDesignBunch(b, None, None, ind_start, ind_stop)
        e_kin_out = b.getSyncParticle().kinEnergy()
        if e_kin_max < e_kin_out:
            e_kin_max = e_kin_out
            phase_max = phase_deg
    cav_phase = phase_max + synchPhaseDeg
    rf_cav.setPhase(cav_phase * math.pi / 180.0)
    # print "debug cav=",cav_name," phase=",cav_phase," delta_e max=",(e_kin_max-e_kin_in)/1.0e-3
    return cav_phase


# -------------------------------------------------------------------
#          Start of the script
# -------------------------------------------------------------------

random.seed(100)

names = ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4", "DTL5", "DTL6", "CCL1", "CCL2", "CCL3", "CCL4", "SCLMed", "SCLHigh", "HEBT1", "HEBT2"]

# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.02)

# ---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_sts_linac.xml"

# ---- make lattice from XML file
accLattice = sns_linac_factory.getLinacAccLattice(names, xml_file_name)

print("Linac lattice is ready. L=", accLattice.getLength())

# ----set up RF Gap Model -------------
# ---- There are three available models at this moment
# ---- BaseRfGap  uses only E0TL*cos(phi)*J0(kr) with E0TL = const
# ---- MatrixRfGap uses a matrix approach like envelope codes
# ---- RfGapTTF uses Transit Time Factors (TTF) like PARMILA
# cppGapModel = BaseRfGap
# cppGapModel = MatrixRfGap
cppGapModel = RfGapTTF
rf_gaps = accLattice.getRF_Gaps()
for rf_gap in rf_gaps:
    rf_gap.setCppGapModel(cppGapModel())

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
# set FFT 3D Space Charge
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
aprtNodes = Add_rfgap_apertures_to_lattice(accLattice, aprtNodes)
aprtNodes = AddMEBTChopperPlatesAperturesToSNS_Lattice(accLattice, aprtNodes)

x_size = 0.042
y_size = 0.042
aprtNodes = AddScrapersAperturesToLattice(accLattice, "MEBT_Diag:H_SCRP", x_size, y_size, aprtNodes)

x_size = 0.042
y_size = 0.042
aprtNodes = AddScrapersAperturesToLattice(accLattice, "MEBT_Diag:V_SCRP", x_size, y_size, aprtNodes)

"""
for node in aprtNodes:
	print "aprt=",node.getName()," pos =",node.getPosition()
"""

# ------ definitions of the activated cavities with Amp = 1
sclHigh_accSeq = accLattice.getSequence("SCLHigh")
cav_names = []
cav_names = [
    "SCL:Cav32a",
    "SCL:Cav32b",
    "SCL:Cav32c",
    "SCL:Cav32d",
] + cav_names
cav_names = [
    "SCL:Cav31a",
    "SCL:Cav31b",
    "SCL:Cav31c",
    "SCL:Cav31d",
] + cav_names
cav_names = [
    "SCL:Cav30a",
    "SCL:Cav30b",
    "SCL:Cav30c",
    "SCL:Cav30d",
] + cav_names
cav_names = [
    "SCL:Cav29a",
    "SCL:Cav29b",
    "SCL:Cav29c",
    "SCL:Cav29d",
] + cav_names
cav_names = [
    "SCL:Cav28a",
    "SCL:Cav28b",
    "SCL:Cav28c",
    "SCL:Cav28d",
] + cav_names
cav_names = [
    "SCL:Cav27a",
    "SCL:Cav27b",
    "SCL:Cav27c",
    "SCL:Cav27d",
] + cav_names
cav_names = [
    "SCL:Cav25a",
    "SCL:Cav25b",
    "SCL:Cav25c",
    "SCL:Cav25d",
] + cav_names

rf_cavs = []
for cav_name in cav_names:
    rf_cav = sclHigh_accSeq.getRF_Cavity(cav_name)
    rf_cav.setAmp(1.0)
    rf_cavs.append(rf_cav)


# ------ definitions synchronous phases for the activated cavities with Amp = 1
cav_synch_phases = {}

"""
cav_synch_phases["SCL:Cav25a"] = -20.0
cav_synch_phases["SCL:Cav25b"] = -20.0
cav_synch_phases["SCL:Cav25c"] = -20.0
cav_synch_phases["SCL:Cav25d"] = -20.0

cav_synch_phases["SCL:Cav27a"] = -30.0
cav_synch_phases["SCL:Cav27b"] = +30.0
cav_synch_phases["SCL:Cav27c"] = -30.0
cav_synch_phases["SCL:Cav27d"] = -0.0

cav_synch_phases["SCL:Cav28a"] = -30.0
cav_synch_phases["SCL:Cav28b"] = +30.0
cav_synch_phases["SCL:Cav28c"] = -30.0
cav_synch_phases["SCL:Cav28d"] = -0.0

cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = +30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -0.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = +30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -0.0

#------scl32 added -----start----
cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = +30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] =  +0.0
#------scl32 added -----stop----


#------scl31 added -----start----
cav_synch_phases["SCL:Cav31a"] = -30.0
cav_synch_phases["SCL:Cav31b"] = +30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = +30.0

cav_synch_phases["SCL:Cav32a"] = -35.0
cav_synch_phases["SCL:Cav32b"] = +35.0
cav_synch_phases["SCL:Cav32c"] = -35.0
cav_synch_phases["SCL:Cav32d"] = -20.0
#------scl31 added -----stop----

#------scl30 added -----start----
cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -30.0

cav_synch_phases["SCL:Cav31a"] = +30.0
cav_synch_phases["SCL:Cav31b"] = -30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = +30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = +30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
#------scl30 added -----stop----


#------scl29 added -----start----
cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = -30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -30.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -10.0
cav_synch_phases["SCL:Cav30d"] = +10.0

cav_synch_phases["SCL:Cav31a"] = -0.0
cav_synch_phases["SCL:Cav31b"] = -0.0
cav_synch_phases["SCL:Cav31c"] = -0.0
cav_synch_phases["SCL:Cav31d"] = -30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
#------scl29 added -----stop----


#------scl28 added -----start----
cav_synch_phases["SCL:Cav28a"] = -30.0
cav_synch_phases["SCL:Cav28b"] = -30.0
cav_synch_phases["SCL:Cav28c"] = -30.0
cav_synch_phases["SCL:Cav28d"] = -30.0

cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = -30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -30.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -30.0

cav_synch_phases["SCL:Cav31a"] = -30.0
cav_synch_phases["SCL:Cav31b"] = -30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = -30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
#------scl28 added -----stop----

#------scl27 added -----start----
cav_synch_phases["SCL:Cav27a"] = -30.0
cav_synch_phases["SCL:Cav27b"] = -30.0
cav_synch_phases["SCL:Cav27c"] = -30.0
cav_synch_phases["SCL:Cav27d"] = -30.0

cav_synch_phases["SCL:Cav28a"] = -30.0
cav_synch_phases["SCL:Cav28b"] = -30.0
cav_synch_phases["SCL:Cav28c"] = -30.0
cav_synch_phases["SCL:Cav28d"] = -30.0

cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = -30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -30.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -30.0

cav_synch_phases["SCL:Cav31a"] = -30.0
cav_synch_phases["SCL:Cav31b"] = -30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = -30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
#------scl27 added -----stop----
"""

# ------scl25 added -----start----
cav_synch_phases["SCL:Cav25a"] = -30.0
cav_synch_phases["SCL:Cav25b"] = -30.0
cav_synch_phases["SCL:Cav25c"] = -30.0
cav_synch_phases["SCL:Cav25d"] = -30.0

cav_synch_phases["SCL:Cav27a"] = -30.0
cav_synch_phases["SCL:Cav27b"] = -30.0
cav_synch_phases["SCL:Cav27c"] = -30.0
cav_synch_phases["SCL:Cav27d"] = -30.0

cav_synch_phases["SCL:Cav28a"] = -30.0
cav_synch_phases["SCL:Cav28b"] = -30.0
cav_synch_phases["SCL:Cav28c"] = -30.0
cav_synch_phases["SCL:Cav28d"] = -30.0

cav_synch_phases["SCL:Cav29a"] = -30.0
cav_synch_phases["SCL:Cav29b"] = -30.0
cav_synch_phases["SCL:Cav29c"] = -30.0
cav_synch_phases["SCL:Cav29d"] = -30.0

cav_synch_phases["SCL:Cav30a"] = -30.0
cav_synch_phases["SCL:Cav30b"] = -30.0
cav_synch_phases["SCL:Cav30c"] = -30.0
cav_synch_phases["SCL:Cav30d"] = -30.0

cav_synch_phases["SCL:Cav31a"] = -30.0
cav_synch_phases["SCL:Cav31b"] = -30.0
cav_synch_phases["SCL:Cav31c"] = -30.0
cav_synch_phases["SCL:Cav31d"] = -30.0

cav_synch_phases["SCL:Cav32a"] = -30.0
cav_synch_phases["SCL:Cav32b"] = -30.0
cav_synch_phases["SCL:Cav32c"] = -30.0
cav_synch_phases["SCL:Cav32d"] = -30.0
# ------scl25 added -----stop----


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
bunch_gen.setBeamCurrent(50.0)

bunch_in = bunch_gen.getBunch(nParticles=100000, distributorClass=WaterBagDist3D)
# bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)
# bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = KVDist3D)

print("Bunch Generation completed.")

for cav_name in cav_names:
    cav_synch_phase = cav_synch_phases[cav_name]
    setSynchPhase(bunch_in, accLattice, cav_name, cav_synch_phase)

# set up design
accLattice.trackDesignBunch(bunch_in)

# ----- Print out the average phases for the RF gaps in SCL RF Cavities
# ----- the charge of H- is negative, so the phases are shifted by -180.
cavs = accLattice.getRF_Cavities()
for cav in cavs:
    avg_phase = phaseNearTargetPhaseDeg(cav.getAvgGapPhaseDeg() - 180.0, 0.0)
    print("debug cav=", cav.getName(), " gaps phaseAvg=", avg_phase, " cav_amp=", cav.getAmp())

print("Design tracking completed.")

# track through the lattice
paramsDict = {"old_pos": -1.0, "count": 0, "pos_step": 0.1}
actionContainer = AccActionsContainer("Test Design Bunch Tracking")


pos_start = 0.0

twiss_analysis = BunchTwissAnalysis()

file_out = open("pyorbit_sts_twiss_sizes_ekin.dat", "w")

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

# ---- If we need to repeat the tracking we can do it starting after RF SCL23d cavity
# ---- Here we do not need this, but the necessary actions are there

# ---- Let's find the last RF gap of SCL23d index in the lattice
rf_cav = accLattice.getRF_Cavity("SCL:Cav23d")
rf_gaps = rf_cav.getRF_GapNodes()
ind_stop = accLattice.getNodeIndex(rf_gaps[len(rf_gaps) - 1])

# ------------Track bunch to the last RF gap in SCL:Cav23d
accLattice.trackBunch(bunch_in, paramsDict=paramsDict, actionContainer=actionContainer, index_start=-1, index_stop=ind_stop)

# ---- Dump and read the bunch into the disk if necessary
# bunch_in.dumpBunch("bunch_sts_after_scl_23d_rg06.dat")

# bunch_in.deleteAllParticles()
# bunch_in.readBunch("bunch_sts_after_scl_23d_rg06.dat")

# ------------Now track the bunch to the end
accLattice.trackBunch(bunch_in, paramsDict=paramsDict, actionContainer=actionContainer, index_start=(ind_stop + 1))

time_exec = time.process_time() - time_start
print("time[sec]=", time_exec)

file_out.close()

"""
#-------------------------------------------------
# Plot the phi-dE phase space at the end
#-------------------------------------------------

z_deg_arr = []
dE_arr = []
for ind in range(bunch_in.getSize()):
	(z_deg,dE) = (bunch_in.z(ind)*bunch_gen.getZtoPhaseCoeff(bunch_in),bunch_in.dE(ind)*1000.)
	z_deg_arr.append(z_deg)
	dE_arr.append(dE)

import Gnuplot
data = Gnuplot.Data(z_deg_arr,dE_arr,with_='p 19', title='Z-dE')
gp = Gnuplot.Gnuplot(persist = 1)
gp('set grid')
gp('set key left')
gp.title("Long. Phase Space")
gp('set xlabel "phi,deg"')
gp('set ylabel "dE, MeV"')
#gp('set pointsize 1.5')
gp.plot(data)
raw_input('Please press return to stop:\n')

sys.exit(1)
"""
