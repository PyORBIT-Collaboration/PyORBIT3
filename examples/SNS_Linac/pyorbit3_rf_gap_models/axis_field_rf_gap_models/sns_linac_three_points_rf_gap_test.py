#! /usr/bin/env python

"""
This script will check the three-point TTF gap model in the SNS linac lattice.
Here we will use only MEBT cavities and DTL1 rf gap 1 (1 and 2),
because we have to specify the exact input energy.
It also includes a speed test and a memory leak test.

Tests include the usual RF Gap class instance extracted from the standard lattice:
rf_gap =  accLattice.getNodeForName("name of the RF gap")

and the new instance of the AxisFieldRF_Gap instance
that uses the RfGapThreePointTTF class instance from C++ linac package:
three_point_gap = AxisFieldRF_Gap(rf_gap)

All C++ tracker classes are from "linac" package:
>from linac import BaseRfGap, MatrixRfGap, RfGapTTF, RfGapThreePointTTF


The script tracks the Bunch with only a synchronous particle at the right energy.
We are going to compare energy gains for different models for MEBT re-bunchers
and the first gap of DTL1.
We also going to print the synch. particle phases at different parts of the
distributed RF field in the AxisFieldRF_Gap and the average phase which should
be almost equal to the synch. particle phase in the zero length BaseRfGap
instance.

"""

import sys
import math
import random
import time
import orbit.core

import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

# from linac import the C++ RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.py_linac.lattice import Drift, AxisFieldRF_Gap

from bunch import Bunch


# names = ["MEBT","DTL1","DTL2","DTL3","DTL4","DTL5","DTL6","CCL1","CCL2","CCL3","CCL4","SCLMed","SCLHigh","HEBT1","HEBT2"]
names = ["MEBT", "DTL1"]

# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.01)

# ---- the XML file name with the structure
xml_file_name = "../../sns_linac_xml/sns_linac.xml"

# ---- make lattice from XML file
accLattice = sns_linac_factory.getLinacAccLattice(names, xml_file_name)

print("Linac lattice is ready. L=", accLattice.getLength())

bunch_init = Bunch()
# set H- mass
# bunch_init.mass(0.9382723 + 2*0.000511)
bunch_init.mass(0.939294)
bunch_init.charge(-1.0)
bunch_init.getSyncParticle().kinEnergy(0.0025)

# ----------------------------------------
# Test MEBT Rb4
# ----------------------------------------

bunch = Bunch()
bunch_init.copyEmptyBunchTo(bunch)
bunch.getSyncParticle().time(0.0)
eKin_init = bunch.getSyncParticle().kinEnergy()

rf_gap = accLattice.getNodeForName("MEBT_RF:Bnch04:Rg01")
three_point_gap = AxisFieldRF_Gap(rf_gap)
print("debug =====================================================================")
print("TTF        gap =", rf_gap.getName())
print("ThreePoint gap =", three_point_gap.getName())


z_step = 0.01
print("debug z_step[mm]=", z_step * 1000.0)

dir_location = "../../sns_rf_fields/"
three_point_gap.readAxisFieldFile(dir_location, three_point_gap.getParam("EzFile"), z_step)

# -------------------------- Track Bunch (Design)---------------------
three_point_gap.trackDesignBunch(bunch)

print("debug ====================Design===========================================")
print("debug gap=", three_point_gap.getName())
print("debug AxisFieldRF_Gap delta_eKin[MeV]=", (bunch.getSyncParticle().kinEnergy() - eKin_init) * 1000.0)
print("debug TTF Gap         delta_eKin[MeV]=", math.fabs(rf_gap.getParam("E0TL") * math.cos(rf_gap.getParam("gap_phase"))) * 1000.0)
print("debug gap_phase=", three_point_gap.getGapPhase() * 180.0 / math.pi)

phase_avg = 0.0
for [pos, phase] in three_point_gap.gap_phase_vs_z_arr:
    print("pos = %6.1f " % (pos * 1000.0), " phase = %6.2f " % (phase * 180.0 / math.pi))
    phase_avg += phase
phase_avg /= len(three_point_gap.gap_phase_vs_z_arr)
print("debug avg_phase = %6.2f " % (phase_avg * 180.0 / math.pi))


# -------------------------- Track Bunch (not Design)-----------------
# we have to create bunch again because trackDesignBunch(bunch) method for node
# does not create an empty copy of the bunch like the lattice method.
# --------------------------------------------------------------------
bunch = Bunch()
bunch_init.copyEmptyBunchTo(bunch)
bunch.getSyncParticle().time(0.0)
eKin_init = bunch.getSyncParticle().kinEnergy()

cav_phase = three_point_gap.getRF_Cavity().getPhase()
cav_phase_shift_deg = -0.0
three_point_gap.getRF_Cavity().setPhase(cav_phase + cav_phase_shift_deg * math.pi / 180.0)


three_point_gap.trackBunch(bunch)
print("debug ===================Not Design========================================")
print("debug cav phase shift = ", cav_phase_shift_deg)
print("debug gap=", three_point_gap.getName())
print("debug  AxisFieldRF_Gap delta_eKin[MeV]=", (bunch.getSyncParticle().kinEnergy() - eKin_init) * 1000.0)
print("debug TTF Gap          delta_eKin[MeV]=", math.fabs(rf_gap.getParam("E0TL") * math.cos(rf_gap.getParam("gap_phase"))) * 1000.0)
print("debug gap_phase=", three_point_gap.getGapPhase() * 180.0 / math.pi)

phase_avg = 0.0
for [pos, phase] in three_point_gap.gap_phase_vs_z_arr:
    print("pos = %6.1f " % (pos * 1000.0), " phase = %6.2f " % (phase * 180.0 / math.pi))
    phase_avg += phase
phase_avg /= len(three_point_gap.gap_phase_vs_z_arr)
print("debug avg_phase = %6.2f " % (phase_avg * 180.0 / math.pi))

# sys.exit(1)

print("#----------------------------------------------------------")
print("#          DTL gaps 1 and 2")
print("#----------------------------------------------------------")

rf_gap1 = accLattice.getNodeForName("DTL_RF:Cav01:Rg01")
rf_gap2 = accLattice.getNodeForName("DTL_RF:Cav01:Rg02")

three_point_gap1 = AxisFieldRF_Gap(rf_gap1)
three_point_gap2 = AxisFieldRF_Gap(rf_gap2)

print("TTF        gap 1 =", rf_gap1.getName())
print("ThreePoint gap 1 =", three_point_gap1.getName())
print("TTF        gap 2 =", rf_gap2.getName())
print("ThreePoint gap 2 =", three_point_gap2.getName())

bunch = Bunch()
bunch_init.copyEmptyBunchTo(bunch)
bunch.getSyncParticle().time(0.0)
eKin_init = bunch.getSyncParticle().kinEnergy()

z_step = 0.01
print("debug z_step[mm]=", z_step * 1000.0)
three_point_gap1.readAxisFieldFile(dir_location, three_point_gap1.getParam("EzFile"), z_step)
three_point_gap2.readAxisFieldFile(dir_location, three_point_gap2.getParam("EzFile"), z_step)

print("debug ========================Design Tracking========================")
three_point_gap1.trackDesignBunch(bunch)
print("debug ===============================================================")
print("debug gap=", three_point_gap1.getName())
print("debug AxisFieldRF_Gap delta_eKin[MeV]=", (bunch.getSyncParticle().kinEnergy() - eKin_init) * 1000.0)
print("debug TTF Gap         delta_eKin[MeV]=", math.fabs(rf_gap1.getParam("E0TL") * math.cos(rf_gap1.getParam("gap_phase"))) * 1000.0)
print("debug gap_phase=", three_point_gap1.getGapPhase() * 180.0 / math.pi)

phase_avg = 0.0
for [pos, phase] in three_point_gap1.gap_phase_vs_z_arr:
    print("pos = %6.1f " % (pos * 1000.0), " phase = %6.2f " % (phase * 180.0 / math.pi))
    phase_avg += phase
phase_avg /= len(three_point_gap1.gap_phase_vs_z_arr)
print("debug avg_phase = %6.2f " % (phase_avg * 180.0 / math.pi))


eKin_init = bunch.getSyncParticle().kinEnergy()
three_point_gap2.trackDesignBunch(bunch)
print("debug ===============================================================")
print("debug gap=", three_point_gap2.getName())
print("debug AxisFieldRF_Gap delta_eKin[MeV]=", (bunch.getSyncParticle().kinEnergy() - eKin_init) * 1000.0)
print("debug TTF Gap         delta_eKin[MeV]=", math.fabs(rf_gap2.getParam("E0TL") * math.cos(rf_gap2.getParam("gap_phase"))) * 1000.0)
print("debug gap_phase=", three_point_gap2.getGapPhase() * 180.0 / math.pi)

phase_avg = 0.0
for [pos, phase] in three_point_gap2.gap_phase_vs_z_arr:
    print("pos = %6.1f " % (pos * 1000.0), " phase = %6.2f " % (phase * 180.0 / math.pi))
    phase_avg += phase
phase_avg /= len(three_point_gap2.gap_phase_vs_z_arr)
print("debug avg_phase = %6.2f " % (phase_avg * 180.0 / math.pi))

print("debug ========================NOT Design Tracking====================")

bunch = Bunch()
bunch_init.copyEmptyBunchTo(bunch)
bunch.getSyncParticle().time(0.0)
eKin_init = bunch.getSyncParticle().kinEnergy()

three_point_gap1.trackBunch(bunch)
print("debug ===============================================================")
print("debug gap=", three_point_gap1.getName())
print("debug AxisFieldRF_Gap delta_eKin[MeV]=", (bunch.getSyncParticle().kinEnergy() - eKin_init) * 1000.0)
print("debug TTF Gap         delta_eKin[MeV]=", math.fabs(rf_gap1.getParam("E0TL") * math.cos(rf_gap1.getParam("gap_phase"))) * 1000.0)
print("debug gap_phase=", three_point_gap1.getGapPhase() * 180.0 / math.pi)

phase_avg = 0.0
for [pos, phase] in three_point_gap1.gap_phase_vs_z_arr:
    print("pos = %6.1f " % (pos * 1000.0), " phase = %6.2f " % (phase * 180.0 / math.pi))
    phase_avg += phase
phase_avg /= len(three_point_gap1.gap_phase_vs_z_arr)
print("debug avg_phase = %6.2f " % (phase_avg * 180.0 / math.pi))


eKin_init = bunch.getSyncParticle().kinEnergy()
three_point_gap2.trackBunch(bunch)
print("debug ===============================================================")
print("debug gap=", three_point_gap2.getName())
print("debug  AxisFieldRF_Gap delta_eKin[MeV]=", (bunch.getSyncParticle().kinEnergy() - eKin_init) * 1000.0)
print("debug TTF Gap          delta_eKin[MeV]=", math.fabs(rf_gap2.getParam("E0TL") * math.cos(rf_gap2.getParam("gap_phase"))) * 1000.0)
print("debug gap_phase=", three_point_gap2.getGapPhase() * 180.0 / math.pi)

phase_avg = 0.0
for [pos, phase] in three_point_gap2.gap_phase_vs_z_arr:
    print("pos = %6.1f " % (pos * 1000.0), " phase = %6.2f " % (phase * 180.0 / math.pi))
    phase_avg += phase
phase_avg /= len(three_point_gap2.gap_phase_vs_z_arr)
print("debug avg_phase = %6.2f " % (phase_avg * 180.0 / math.pi))

sys.exit(0)

# ------------------------------------
#       SPEED and MEMORY LEAK TEST
# ------------------------------------


three_point_gap = three_point_gap1
rf_gap = rf_gap1

paramsDict = {}

count = 0
bunch = Bunch()
time_start = time.time()
eKin_init = bunch_init.getSyncParticle().kinEnergy()

while 1 < 2:
    count += 1
    bunch_init.copyEmptyBunchTo(bunch)
    bunch.getSyncParticle().time(0.0)
    # three_point_gap.trackDesignBunch(bunch)
    three_point_gap.trackBunch(bunch)
    if count % 100 == 0:
        speed = count / (time.time() - time_start)
        print("debug ================  count=", count, " speed [calc/sec] =", speed)
        print("debug         delta_eKin[MeV]=", (bunch.getSyncParticle().kinEnergy() - eKin_init) * 1000.0)
        print("debug TTF Gap delta_eKin[MeV]=", math.fabs(rf_gap.getParam("E0TL") * math.cos(rf_gap.getParam("gap_phase"))) * 1000.0)
