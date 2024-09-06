#! /usr/bin/env python

"""
This script will track the bunch through the SNS MEBT + DTL Linac.
The bunch will be modified at the beginning and the results of tracking
(bunch sizes) should be the same.
"""

import os
import math
import random
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from orbit.core import orbit_mpi
from orbit.core.orbit_mpi import mpi_comm, mpi_datatype

# from linac import the C++ RF gap classes
from orbit.core.linac import BaseRfGap

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D

from orbit.core.bunch import Bunch, BunchTwissAnalysis, SynchPartRedefinitionZdE

from orbit.lattice import AccNode, AccActionsContainer

from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice

from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes
from orbit.py_linac.lattice_modifications import Replace_Quads_to_OverlappingQuads_Nodes

from orbit.py_linac.overlapping_fields import SNS_EngeFunctionFactory

from orbit.py_linac.lattice import LinacPhaseApertureNode


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


class SNS_Linac_BunchGenerator:
    """
    Generates the pyORBIT SNS Linac Bunches.
    Twiss parameters has the fol following units: x in [m], xp in [rad]
    and the X and Y emittances are un-normalized. The longitudinal emittance
    is in [GeV*m].
    """

    def __init__(self, twissX, twissY, twissZ, frequency=402.5e6):
        self.twiss = (twissX, twissY, twissZ)
        self.bunch_frequency = frequency
        self.bunch = Bunch()
        syncPart = self.bunch.getSyncParticle()
        # set H- mass
        # self.bunch.mass(0.9382723 + 2*0.000511)
        self.bunch.mass(0.939294)
        self.bunch.charge(-1.0)
        syncPart.kinEnergy(0.0025)
        self.c = 2.99792458e8  # speed of light in m/sec
        self.beam_current = 38.0  # beam current in mA , design = 38 mA
        self.rf_wave_lenght = self.c / self.bunch_frequency
        self.si_e_charge = 1.6021773e-19

    def getKinEnergy(self):
        """
        Returns the kinetic energy in GeV
        """
        return self.bunch.getSyncParticle().kinEnergy()

    def setKinEnergy(self, e_kin=0.0025):
        """
        Sets the kinetic energy in GeV
        """
        self.bunch.getSyncParticle().kinEnergy(e_kin)

    def getZtoPhaseCoeff(self, bunch):
        """
        Returns the coefficient to calculate phase in degrees from the z-coordinate.
        """
        bunch_lambda = bunch.getSyncParticle().beta() * self.rf_wave_lenght
        phase_coeff = 360.0 / bunch_lambda
        return phase_coeff

    def getBeamCurrent(self):
        """
        Returns the beam currect in mA
        """
        return self.beam_current

    def setBeamCurrent(self, current):
        """
        Sets  the beam currect in mA
        """
        self.beam_current = current

    def getBunch(self, nParticles=0, distributorClass=WaterBagDist3D, cut_off=-1.0):
        """
        Returns the pyORBIT bunch with particular number of particles.
        """
        comm = mpi_comm.MPI_COMM_WORLD
        rank = orbit_mpi.MPI_Comm_rank(comm)
        size = orbit_mpi.MPI_Comm_size(comm)
        data_type = mpi_datatype.MPI_DOUBLE
        main_rank = 0
        bunch = Bunch()
        self.bunch.copyEmptyBunchTo(bunch)
        macrosize = self.beam_current * 1.0e-3 / self.bunch_frequency
        macrosize /= math.fabs(bunch.charge()) * self.si_e_charge
        distributor = None
        if distributorClass == WaterBagDist3D:
            distributor = distributorClass(self.twiss[0], self.twiss[1], self.twiss[2])
        else:
            distributor = distributorClass(self.twiss[0], self.twiss[1], self.twiss[2], cut_off)
        bunch.getSyncParticle().time(0.0)
        for i in range(nParticles):
            (x, xp, y, yp, z, dE) = distributor.getCoordinates()
            (x, xp, y, yp, z, dE) = orbit_mpi.MPI_Bcast((x, xp, y, yp, z, dE), data_type, main_rank, comm)
            if i % size == rank:
                bunch.addParticle(x, xp, y, yp, z, dE)
        nParticlesGlobal = bunch.getSizeGlobal()
        bunch.macroSize(macrosize / nParticlesGlobal)
        return bunch


script_dir = os.path.dirname(__file__)
random.seed(100)

names = ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4", "DTL5", "DTL6"]
names = [
    "MEBT",
]
names = ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4"]


# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

# ---- the XML file name with the structure
xml_file_name = os.path.join(script_dir, "../sns_linac_xml/sns_linac.xml")

# ---- make lattice from XML file
accLattice = sns_linac_factory.getLinacAccLattice(names, xml_file_name)

print("Linac lattice is ready. L=", accLattice.getLength())

# ----set up RF Gap Model -------------
# ---- There are three available models at this moment
# ---- BaseRfGap  uses only E0TL*cos(phi)*J0(kr) with E0TL = const
# ---- MatrixRfGap uses a matrix approach like envelope codes
# ---- RfGapTTF uses Transit Time Factors (TTF) like PARMILA
cppGapModel = BaseRfGap
# cppGapModel = BaseRfGap_slow
# cppGapModel = MatrixRfGap
# cppGapModel = RfGapTTF
rf_gaps = accLattice.getRF_Gaps()
for rf_gap in rf_gaps:
    rf_gap.setCppGapModel(cppGapModel())

# ------------------------------------------------------------------
# ---- BaseRF_Gap and Quads will be replaced for specified sequences
# ------------------------------------------------------------------

# ---- longitudinal step along the distributed fields lattice
z_step = 0.005

# ---- axis fields files location
dir_location = os.path.join(script_dir, "../../../../examples/SNS_Linac/sns_rf_fields/")

# Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,["MEBT",])

Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(
    accLattice, z_step, dir_location, ["MEBT", "DTL1", "DTL2", "DTL3", "DTL4"], [], SNS_EngeFunctionFactory
)
# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["SCLMed","SCLHigh"],[],SNS_EngeFunctionFactory)
Replace_Quads_to_OverlappingQuads_Nodes(
    accLattice,
    z_step,
    [
        "MEBT",
    ],
    [],
    SNS_EngeFunctionFactory,
)
# Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,["MEBT","DTL1"],[],SNS_EngeFunctionFactory)
# Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,["DTL1",],[],SNS_EngeFunctionFactory)

# ---- setup the the linac specific trackers
accLattice.setLinacTracker(True)


frequency = 402.5e6
node_pos_dict = accLattice.getNodePositionsDict()
rf_gaps = accLattice.getRF_Gaps()
aprtNodes = []
for rf_gap in rf_gaps:
    phaseAperture = LinacPhaseApertureNode(frequency, rf_gap.getName() + ":phaseAprt")
    pos = node_pos_dict[rf_gap][0]
    phaseAperture.setPosition(pos)
    phaseAperture.setMinMaxPhase(-180.0 * 2, +180.0 * 2)
    rf_gap.addChildNode(phaseAperture, AccNode.EXIT)
    aprtNodes.append(phaseAperture)


# -----------------------------------------------------
# Set up Space Charge Acc Nodes
# -----------------------------------------------------
from orbit.space_charge.sc3d import setUniformEllipsesSCAccNodes
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse

sc_path_length_min = 0.02

print("Set up Space Charge nodes. ")

# set of uniformly charged ellipses Space Charge
nEllipses = 1
calcUnifEllips = SpaceChargeCalcUnifEllipse(nEllipses)
space_charge_nodes = setUniformEllipsesSCAccNodes(accLattice, sc_path_length_min, calcUnifEllips)

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

bunch_in = bunch_gen.getBunch(nParticles=100000, distributorClass=WaterBagDist3D)
# bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)
# bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = KVDist3D)

print("Bunch Generation completed.")

# set up design
accLattice.trackDesignBunch(bunch_in)

print("Design tracking completed.")

delta_E = 0.0000
bunch_tmp = Bunch()
bunch_in.copyBunchTo(bunch_tmp)

redefCalc = SynchPartRedefinitionZdE()
redefCalc.shift_dE(bunch_tmp, delta_E)

# ---------------------------------------
redefCalc.analyzeBunch(bunch_tmp)
dE_avg = redefCalc.getAvg_dE()
eKin = bunch_tmp.getSyncParticle().kinEnergy()
print("debug =============== shift E done before tracking =========")
print("debug dE_avg=", dE_avg)
print("debug eKin  =", eKin)
print("debug dE_avg + eKin =", (dE_avg + eKin))
print("debug =======================================================")

# ----------------------------------------------
#  Bunch phase shifting
# ----------------------------------------------
bunch_lambda = bunch_in.getSyncParticle().beta() * 2.99792458e8 / 402.5e6
phase_coeff = 360.0 / bunch_lambda
delta_phase = 20.0  # deg
delta_z = delta_phase / phase_coeff
delta_t = (delta_phase / 360.0) * 1.0 / 402.5e6

print("debug delta_phi[deg]=", delta_phase)
print("debug    delta_z[mm]=", delta_z * 1000.0)
redefCalc.shift_Z(bunch_tmp, delta_z)

# -------------------------------------------------------------------
# track through the lattice
paramsDict = {"old_pos": -1.0, "count": 0, "pos_step": 0.01}
actionContainer = AccActionsContainer("Bunch Tracking")

pos_start = 0.0

twiss_analysis = BunchTwissAnalysis()

file_out = open("pyorbit_twiss_sizes_ekin.dat", "w")

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
    redefCalc.analyzeBunch(bunch)
    eKin += redefCalc.getAvg_dE() * 1.0e3
    s = " %50s  %4.5f " % (node.getName(), pos + pos_start)
    s += "   %6.4f  %6.4f  %6.4f  %6.4f   " % (alphaX, betaX, emittX, norm_emittX)
    s += "   %6.4f  %6.4f  %6.4f  %6.4f   " % (alphaY, betaY, emittY, norm_emittY)
    s += "   %6.4f  %6.4f  %6.4f  %6.4f   " % (alphaZ, betaZ, emittZ, phi_de_emittZ)
    s += "   %5.3f  %5.3f  %5.3f " % (x_rms, y_rms, z_rms_deg)
    s += "  %10.6f   %8d " % (eKin, nParts)
    file_out.write(s + "\n")
    file_out.flush()
    s_prt = " %5d  %50s  %4.5f " % (paramsDict["count"], node.getName(), pos + pos_start)
    s_prt += "  %5.3f  %5.3f   %5.3f " % (x_rms, y_rms, z_rms_deg)
    s_prt += "  %10.6f   %8d " % (eKin, nParts)
    print(s_prt)


def action_exit(paramsDict):
    action_entrance(paramsDict)


actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

time_start = time.process_time()

accLattice.trackBunch(bunch_tmp, paramsDict=paramsDict, actionContainer=actionContainer)

time_exec = time.process_time() - time_start
print("time[sec]=", time_exec)

file_out.close()


redefCalc.analyzeBunch(bunch_tmp)
dE_avg = redefCalc.getAvg_dE()
eKin = bunch_tmp.getSyncParticle().kinEnergy()
print("debug dE_avg=", dE_avg)
print("debug eKin  =", eKin)
print("debug dE_avg + eKin =", (dE_avg + eKin))

out = read_lines("pyorbit_twiss_sizes_ekin.dat")
expected_out = read_lines(os.path.join(script_dir, "expected_pyorbit_twiss_sizes_ekin.dat"))


def test_file():
    result = out == expected_out
    assert result


os.remove("pyorbit_twiss_sizes_ekin.dat")
