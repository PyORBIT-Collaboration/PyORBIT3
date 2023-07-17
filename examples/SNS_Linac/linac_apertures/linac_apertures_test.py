#!/usr/bin/env python

# --------------------------------------------------------
# Test of the Linac Aperture for quads, RF gaps, and scrapers.
# The example includes quads with overlapping fields and the lost bunch.
# The user can comment out the parts with these quads or the lost bunch
# in the paramsDict dictionary.
# --------------------------------------------------------

import math
import sys
import os

from orbit.core.orbit_mpi import mpi_comm, mpi_datatype, mpi_op, MPI_Comm_rank, MPI_Comm_size, MPI_Bcast

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_to_AxisField_Nodes
from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes
from orbit.py_linac.lattice_modifications import Replace_Quads_to_OverlappingQuads_Nodes

from orbit.py_linac.overlapping_fields import SNS_EngeFunctionFactory

from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice
from orbit.py_linac.lattice_modifications import Add_rfgap_apertures_to_lattice
from orbit.py_linac.lattice_modifications import GetLostDistributionArr
from orbit.py_linac.lattice_modifications import AddScrapersAperturesToLattice
from orbit.py_linac.lattice_modifications import AddMEBTChopperPlatesAperturesToSNS_Lattice

from orbit.lattice import AccLattice, AccActionsContainer

from orbit.py_linac.lattice import BaseLinacNode

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D
from orbit.bunch_generators import TwissAnalysis

from orbit.core.bunch import Bunch, BunchTwissAnalysis

# ---- These are convinience class and function for
# ---- particles attributes: Ids and Initial Coordinates
# ---- They will assign these values for the initial
# ---- bunch and we should see them in the lost bunch
from orbit.bunch_utils import ParticleIdNumber
from orbit.core.orbit_utils import bunch_utils_functions


# --------------------------------------------------------
# This is an example of the custom AccNode subclass
# that can be added to the lattice or attached as a child
# to any BaseLinacNode instance in the lattice.
# --------------------------------------------------------
class BunchDumpNode(BaseLinacNode):
    """
    The class BunchDumpNode writes the information from the bunch
    to the file.
    """

    def __init__(self, name="BunchDump", file_name="bunch.dat"):
        BaseLinacNode.__init__(self, name)
        self.file_name = file_name

    def setFileName(self, file_name):
        self.file_name = file_name

    def getFileName(self):
        return self.file_name

    def track(self, paramsDict):
        if "bunch" in paramsDict:
            bunch = paramsDict["bunch"]
            bunch.dumpBunch(self.file_name)

    def trackDesign(self, paramsDict):
        """
        This method does nothing for this class.
        """
        pass


# -------------------------------------------
# START of Script
# -------------------------------------------

comm = mpi_comm.MPI_COMM_WORLD
rank = MPI_Comm_rank(comm)
size = MPI_Comm_size(comm)
data_type = mpi_datatype.MPI_DOUBLE
main_rank = 0


names = [
    "MEBT",
]

py_orbit_sns_home = "../"

# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

# ---- the XML file name with the structure
xml_file_name = py_orbit_sns_home + "sns_linac_xml/sns_linac.xml"

# ---- make lattice from XML file
accLattice = sns_linac_factory.getLinacAccLattice(names, xml_file_name)

print("1. Linac lattice is ready. L=", accLattice.getLength())

# --------------------------------------------------------
# Replace overlapping quads and/or RF gaps with the special nodes
# --------------------------------------------------------
# ---- z_step defines the tracking (integration) longitudinal step over the
# ---- quads' and RF gaps' fields
z_step = 0.001

# ---- RF gaps' axis field files location
dir_location = "../sns_rf_fields/"

# Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,["MEBT",])

# Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes(accLattice,z_step,dir_location,["MEBT",],[],SNS_EngeFunctionFactory)

Replace_Quads_to_OverlappingQuads_Nodes(
    accLattice,
    z_step,
    [
        "MEBT",
    ],
    [],
    SNS_EngeFunctionFactory,
)

print("2. Linac lattice is ready. L=", accLattice.getLength())

node_pos_dict = accLattice.getNodePositionsDict()
nodes = accLattice.getNodes()
for node in nodes:
    (posBefore, posAfter) = node_pos_dict[node]
    print("%45s      (start stop) L = (%10.4f %10.4f)    %10.4f " % (node.getName(), posBefore, posAfter, (posAfter - posBefore)))

# ----------------------------------------------------------
#      1. add Aperture nodes to the quads in the linac lattice
#      2. add Aperture nodes to the RF gaps in the linac lattice
#      3. add Aperture nodes to the MEBT chopper entrance/exit plates
#         in the linac lattice
#      4. add Aperture nodes to the MEBT scrapers (H and V)
#      5. add a BunchDumpNode to one of the chopper plate
# ----------------------------------------------------------

print("===== Aperture Nodes =======")

aprtNodes = Add_quad_apertures_to_lattice(accLattice)

aprtNodes = Add_rfgap_apertures_to_lattice(accLattice, aprtNodes)

aprtNodes = AddMEBTChopperPlatesAperturesToSNS_Lattice(accLattice, aprtNodes)

x_size = 0.042
y_size = 0.042
aprtNodes = AddScrapersAperturesToLattice(accLattice, "MEBT_Diag:H_SCRP", x_size, y_size, aprtNodes)

x_size = 0.042
y_size = 0.042
aprtNodes = AddScrapersAperturesToLattice(accLattice, "MEBT_Diag:V_SCRP", x_size, y_size, aprtNodes)

# ---- print all aperture nodes and their positions
for node in aprtNodes:
    print("aprt=", node.getName(), " pos =", node.getPosition())

# ---- set the Bunch Dump Node to see the one of the apertures effect
chopper_plate_entr_node = accLattice.getNodeForName("MEBT:ChpPlt:Entr")
aprt_node_1 = None
for node in aprtNodes:
    if node.getName().find("MEBT:ChpPlt:Entr") >= 0:
        aprt_node_1 = node
        break
if chopper_plate_entr_node != None:
    dumpNode = BunchDumpNode()
    dumpNode.setFileName("bunch_at_chpplt.dat")
    dumpNode.setSequence(chopper_plate_entr_node.getSequence())
    chopper_plate_entr_node.addChildNode(dumpNode, node.EXIT)

# -----------------------------------------------------------
#    Bunch Generation
# -----------------------------------------------------------
bunch = Bunch()
# set H- mass
bunch.mass(0.9382723 + 2 * 0.000511)
bunch.charge(-1.0)
bunch.getSyncParticle().kinEnergy(0.0025)

# electron charge in SI
si_e_charge = 1.6021773e-19
frequency = 402.5e6
beam_current = 0.038

N_particles = 10000
macrosize = beam_current / frequency
macrosize /= math.fabs(bunch.charge()) * si_e_charge
macrosize /= N_particles

(alphaX, betaX, emittX) = (-1.39, 0.126, 3.67 * 1.0e-6)
(alphaY, betaY, emittY) = (2.92, 0.281, 3.74 * 1.0e-6)
(alphaZ, betaZ, emittZ) = (0.0, 117.0, 0.0166 * 1.0e-6)

# ---- we artificially increase the emittances to see apertures effects
emittX *= 5.0
emittY *= 10.0

twissX = TwissContainer(alphaX, betaX, emittX)
twissY = TwissContainer(alphaY, betaY, emittY)
twissZ = TwissContainer(alphaZ, betaZ, emittZ)

distributor = WaterBagDist3D(twissX, twissY, twissZ)

for ind in range(N_particles):
    (x, xp, y, yp, z, dE) = distributor.getCoordinates()
    (x, xp, y, yp, z, dE) = MPI_Bcast((x, xp, y, yp, z, dE), data_type, main_rank, comm)
    if ind % size == rank:
        bunch.addParticle(x, xp, y, yp, z, dE)

nParticlesGlobal = bunch.getSizeGlobal()
if rank == 0:
    print("total number of particles =", nParticlesGlobal)
bunch.macroSize(macrosize)

# ----- Assign the Ids to each particle
ParticleIdNumber.addParticleIdNumbers(bunch)
# ---- Copy the coordinates to the initial coordinates
# ---- array of Particles Attribute
bunch_utils_functions.copyCoordsToInitCoordsAttr(bunch)

# set up design
accLattice.trackDesignBunch(bunch)

paramsDict = {}
lost_parts_bunch = Bunch()
paramsDict["lostbunch"] = lost_parts_bunch

actionContainer = AccActionsContainer("Test Apertures")

twiss_analysis = BunchTwissAnalysis()


def action_entrance(paramsDict):
    if isinstance(paramsDict["parentNode"], AccLattice):
        node = paramsDict["node"]
        pos = paramsDict["path_length"]
        bunch = paramsDict["bunch"]
        twiss_analysis.analyzeBunch(bunch)
        x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1] * twiss_analysis.getTwiss(0)[3]) * 1000.0
        y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1] * twiss_analysis.getTwiss(1)[3]) * 1000.0
        z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1] * twiss_analysis.getTwiss(2)[3]) * 1000.0
        eKin = bunch.getSyncParticle().kinEnergy() * 1.0e3
        s = " %35s  %4.5f    %5.3f  %5.3f   %5.3f     %10.6f   %8d " % (
            node.getName(),
            pos,
            x_rms,
            y_rms,
            z_rms,
            eKin,
            bunch.getSizeGlobal(),
        )
        print(s)


actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
accLattice.trackBunch(bunch, paramsDict=paramsDict, actionContainer=actionContainer)

lost_parts_bunch.dumpBunch("lostbunch.dat")

# ----------------------------------------------------------------------------
# This is an example how to analyze the lost bunch with the special function
# GetLostDistributionArr(...)
# ----------------------------------------------------------------------------

aprtNodes_loss_arr = GetLostDistributionArr(aprtNodes, lost_parts_bunch)
total_loss = 0.0
for [aprtNode, loss] in aprtNodes_loss_arr:
    print("aprt. node= %30s " % aprtNode.getName(), " pos= %9.3f " % aprtNode.getPosition(), " loss= %6.0f " % loss)
    total_loss += loss
print("Total loss=", total_loss)
