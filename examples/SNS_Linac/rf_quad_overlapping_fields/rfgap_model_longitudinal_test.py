#--------------------------------------------------------------------
# This script tests the functionality of AxisFieldRF_Gap class.
# It also performs benchmark this calss with the zero-length BaseRfGap() model.
# There are cases:
# =============================================================
# Case 1: Cavity phase scan assuming the phase at z=0 position.
# Case 2: Cavity phase scan assuming the phase at the entrance 
#  of the cavity
# Case 3: Cavity phase scan with non-empty bunch assuming the phase 
#  at the entrance of the cavity
# Case 4: Thin (zero length) RF Gap and Axis Field Gap comparison
#---------------------------------------------------------------------

import math
import random
import time
import sys

# from linac import the C++ RF gap classes
from orbit.core.linac import BaseRfGap, MatrixRfGap, RfGapTTF

from orbit.core.bunch import Bunch
from orbit.py_linac.lattice import BaseRF_Gap, AxisFieldRF_Gap, RF_Cavity
from orbit.py_linac.lattice import Drift

from orbit.core.orbit_utils import Function, SplineCH, GaussLegendreIntegrator, Polynomial
from orbit.utils.fitting import PolynomialFit

from orbit.utils import phaseNearTargetPhase
from orbit.utils import phaseNearTargetPhaseDeg

#--------------------------------------------
#  Auxiliary Functions
#--------------------------------------------

def getRectangularRfAxisField(z_min,z_max):
    """
    The RF axial field step function normalized to the integral value of 1.
    The center position will be assumed at z=0.
    """
    value = 1./(z_max - z_min)
    npoints = 10
    z_step = (z_max - z_min)/(npoints - 1)
    func = Function()
    for ind in range(npoints):
        z = z_min + z_step*ind
        func.add(z,value)
    return func
    
def transitTimeFactor(bunch,cav_frequency,length):
    """
    Returns the TTF parameter for the thin RF Gap 
    with a rectangular axis field.
    beta - v/c -relativistic factor
    w = cav_frequency [Hz]
    L2 = length/2 [m]
    lambda = beta*c/w
    TTF = lambda/L2*sin((L2/lambda)
    """
    L2 = length/2
    beta = bunch.getSyncParticle().beta()
    w = 2*math.pi*cav_frequency
    lmbd = beta*2.997924e+8/w
    TTF = (lmbd/L2)*math.sin(L2/lmbd)
    return TTF
#--------------------------------------------
#  Script start
#--------------------------------------------

#---- The RF axial field step function normalized to the integral value of 1.
z_min = -0.1
z_max = 0.1
axis_field_func = getRectangularRfAxisField(z_min,z_max)

#---- Integral E*L= 20 MeV , here E0L parameter is in GeV
E0L = 0.020
#---- The RF gap object representing one RF gap
accNode = BaseRF_Gap("BaseRfGap")
three_point_gap = AxisFieldRF_Gap(accNode)
three_point_gap.setAsFirstRFGap(True)
three_point_gap.setAxisFieldFunction(axis_field_func, z_step = 0.01)
three_point_gap.setParam("E0L", E0L)
#---- mode = 0 - phase as it is, mode = 1 shift phase by +PI
three_point_gap.addParam("mode", 1)

#---- The RF cavity. Here it has only one RF gap
cav_amp = 1.0
cav_frequency = 704.42e6
cav = RF_Cavity("AxisFieldCavity")
cav.setAmp(cav_amp)
cav.setFrequency(cav_frequency)
cav.setPosition((z_min+z_max)/2.)
cav.addRF_GapNode(three_point_gap)

#---- At the entrance of ESS spoke.
#---- Kinetic energy and mass in GeV
mass = 0.93827208943
Ekin_init = 0.6201
bunch_init = Bunch()
bunch_init.mass(mass)
bunch_init.charge(1.0)
bunch_init.getSyncParticle().kinEnergy(Ekin_init)

#---- Parameters of the cavity phase scans.
#---- Phases are in radians.
Nph = 72
phase_step = 2*math.pi/Nph

#--------------------------------------------------------------------
#  Case 1: Cavity phase scan assuming the phase at z=0 position.
#  We have several types of parameters:
#  cav_phase - the input parameter. It is a phase at 
#              the center of the gap
#  cav_1st_gap_entr_phase0 - the phase at the entrance of the 1-st gap
#              (we have only one gap) for the cav_phase = 0.
#  cav_1st_gap_entr_phase  - the phase at the entrance of the 1-st gap
#              for the defined cav_phase
#  gap_phase - the calculated phase in the center of the RF gap after fitting
#              process to get the cav_phase (default accuracy is 0.001 deg)
#  delta_phase = cav_1st_gap_entr_phase - cav_phase - cav_1st_gap_entr_phase0
#              this phase will be 0 at cav_phase = 0, and it characterizes 
#              how wrong we are from the realistic phase scan when we
#              we define cavity phase at the entrance of the cavity
#  delta_eKin_out - the bunch energy gain after acceleration
#--------------------------------------------------------------------
three_point_gap.getRF_Cavity().setPhase(0.)
bunch = Bunch()
bunch_init.copyEmptyBunchTo(bunch)
three_point_gap.trackDesignBunch(bunch)
cav_1st_gap_entr_phase0 = three_point_gap.getRF_Cavity().getFirstGapEtnrancePhase()

print ("=============================================================")
print ("Case 1: Cavity phase scan assuming the phase at z=0 position.")
print ("=============================================================")

st = "CavPhase[deg]  CavPhaseErr[deg]  DeltaPhase[deg] 1stGapPhase[deg] DeltaE[MeV] "
print (st)

for ind in range(Nph+1):
    cav_phase = ind*phase_step
    three_point_gap.getRF_Cavity().setPhase(cav_phase)
    
    bunch = Bunch()
    bunch_init.copyEmptyBunchTo(bunch)
    three_point_gap.trackDesignBunch(bunch)
    
    cav_1st_gap_entr_phase = phaseNearTargetPhase(three_point_gap.getRF_Cavity().getFirstGapEtnrancePhase(),0.)
    delta_phase = phaseNearTargetPhase(cav_1st_gap_entr_phase - cav_phase - cav_1st_gap_entr_phase0,0.)
    gap_phase = phaseNearTargetPhase(three_point_gap.getGapPhase(),cav_phase)
    delta_eKin_out = bunch.getSyncParticle().kinEnergy() - Ekin_init
    st  = " %+8.2f "%(cav_phase*180./math.pi)
    st += " %+9.6f "%((gap_phase-cav_phase)*180./math.pi)
    st += " %+8.4f "%(delta_phase*180./math.pi)
    st += " %+8.2f "%(cav_1st_gap_entr_phase*180./math.pi)
    st += " %+10.6f "%(delta_eKin_out*1000.)
    print (st)


#--------------------------------------------------------------------
#  Case 2: Cavity phase scan assuming the phase at the entrance 
#  of the cavity
#  We have several types of parameters:
#  cav_phase - the input parameter. It is a phase at 
#              the center of the gap
#  cav_1st_gap_entr_phase0 - the phase at the entrance of the 1-st gap
#              (we have only one gap) for the cav_phase = 0.
#  cav_1st_gap_entr_phase  - the phase at the entrance of the 1-st gap
#              for the defined cav_phase
#  gap_phase - the calculated phase in the center of the RF gap after fitting
#              process to get the cav_phase (default accuracy is 0.001 deg)
#  delta_phase = cav_1st_gap_entr_phase - cav_phase - cav_1st_gap_entr_phase0
#              this phase will be 0 at cav_phase = 0, and it characterizes 
#              how wrong we are from the realistic phase scan when we
#              we define cavity phase at the entrance of the cavity
#  delta_eKin_out - the bunch energy gain after acceleration
#--------------------------------------------------------------------

#---- Set the cavity property
three_point_gap.getRF_Cavity().setUsePhaseAtEntrance(True)

three_point_gap.getRF_Cavity().setPhase(0.)
bunch = Bunch()
bunch_init.copyEmptyBunchTo(bunch)
three_point_gap.trackDesignBunch(bunch)
cav_1st_gap_entr_phase0 = three_point_gap.getRF_Cavity().getFirstGapEtnrancePhase()

print ("=============================================================")
print ("Case 2: Cavity phase scan assuming the phase at the entrance.")
print ("=============================================================")

st = "CavPhase[deg]  CavPhaseErr[deg]  DeltaPhase[deg] 1stGapPhase[deg] DeltaE[MeV] "
print (st)

for ind in range(Nph+1):
    cav_phase = ind*phase_step
    three_point_gap.getRF_Cavity().setPhase(cav_phase)
    
    bunch = Bunch()
    bunch_init.copyEmptyBunchTo(bunch)
    three_point_gap.trackDesignBunch(bunch)
    
    cav_1st_gap_entr_phase = phaseNearTargetPhase(three_point_gap.getRF_Cavity().getFirstGapEtnrancePhase(),0.)
    delta_phase = phaseNearTargetPhase(cav_1st_gap_entr_phase - cav_phase - cav_1st_gap_entr_phase0,0.)
    gap_phase = phaseNearTargetPhase(three_point_gap.getGapPhase(),cav_phase)
    delta_eKin_out = bunch.getSyncParticle().kinEnergy() - Ekin_init
    st  = " %+8.2f "%(cav_phase*180./math.pi)
    st += " %+8.4f "%((gap_phase-cav_phase)*180./math.pi)
    st += " %+8.2f "%(delta_phase*180./math.pi)
    st += " %+8.2f "%(cav_1st_gap_entr_phase*180./math.pi)
    st += " %+10.6f "%(delta_eKin_out*1000.)
    print (st)

#--------------------------------------------------------------------
#  Case 3: Cavity phase scan with non-empty bunch assuming the phase 
#  at the entrance of the cavity.
#  We have several types of parameters:
#  cav_phase - the input parameter. It is a phase at 
#              the center of the gap
#  delta_eKin_out - the bunch energy gain after acceleration
#  (x,xp,y,yp,z,dE) - 6D coordinates of the particle
#--------------------------------------------------------------------

#---- Set the cavity property
three_point_gap.getRF_Cavity().setUsePhaseAtEntrance(True)

three_point_gap.getRF_Cavity().setPhase(0.)
bunch = Bunch()
bunch_init.copyEmptyBunchTo(bunch)
three_point_gap.trackDesignBunch(bunch)

print ("=============================================================")
print ("Case 3: Cavity phase scan with one particle in the bunch.")
print ("=============================================================")

st = "CavPhase[deg]  DeltaE[MeV]  x[mm]  xp[mrad] y[mm] yp[mrad]  z[mm] dE[MeV] "
print (st)

for ind in range(Nph+1):
    cav_phase = ind*phase_step
    three_point_gap.getRF_Cavity().setPhase(cav_phase)
    
    bunch = Bunch()
    bunch_init.copyBunchTo(bunch)
    bunch.addParticle(0.01,0.,0.01,0.,0.,0.)

    three_point_gap.trackBunch(bunch)

    delta_eKin_out = bunch.getSyncParticle().kinEnergy() - Ekin_init
    st  = " %+8.2f  %+8.4f  "%(cav_phase*180./math.pi,delta_eKin_out*1000.)
    st += "  %8.3f %+8.3f "%(bunch.x(0)*1000.,bunch.xp(0)*1000.)
    st += "  %8.3f %+8.3f "%(bunch.y(0)*1000.,bunch.yp(0)*1000.)
    st += "  %8.3f %+8.3f "%(bunch.z(0)*1000.,bunch.dE(0)*1000.)
    print (st)

print ("====================================================================")
print ("Case 4: Thin (zero length) RF Gap and Axis Field Gap comparison")
print ("====================================================================")
TTF = transitTimeFactor(bunch_init,cav_frequency,z_max-z_min)
print ("Zero-length RF Gap: Transit Time Factor TTF=",TTF)

base_rf_gap = BaseRF_Gap("BaseRfGap")
base_rf_gap.setParam("E0TL", E0L*TTF)
#---- setting RF gap tracker
base_rf_gap.setCppGapModel(BaseRfGap())

#---- The RF cavity for zero-length RF gap. Here it has only one RF gap
cav_amp = 1.0
cav_frequency = 704.42e6
cav_base = RF_Cavity("BaseGapCavity")
cav_base.setAmp(cav_amp)
cav_base.setFrequency(cav_frequency)
cav_base.setPosition((z_min+z_max)/2.)
cav_base.addRF_GapNode(base_rf_gap)

#---- Tracking design for both types of RF gaps
base_rf_gap.getRF_Cavity().setPhase(0.)
bunch = Bunch()
bunch_init.copyEmptyBunchTo(bunch)
base_rf_gap.trackDesignBunch(bunch)

three_point_gap.getRF_Cavity().setPhase(0.)
bunch = Bunch()
bunch_init.copyEmptyBunchTo(bunch)
three_point_gap.trackDesignBunch(bunch)

#---- drifts before and after zero-length RF gap
drift_before = Drift("drift_before")
drift_before.setLength((z_max-z_min)/2)
drift_after = Drift("drift_after")
drift_after.setLength((z_max-z_min)/2)

print ("==== 0 - for zero-length and 1 - for axis field RF gap ====") 
st = "CavPhase[deg]  deltaE_0[MeV]  deltaE_1[MeV]  delta_0-1[MeV] delta_phase_0-1[deg]"
print (st)

for ind in range(Nph+1):
    cav_phase = ind*phase_step
    base_rf_gap.getRF_Cavity().setPhase(cav_phase)
    three_point_gap.getRF_Cavity().setPhase(cav_phase)
    
    bunch = Bunch()
    bunch_init.copyBunchTo(bunch)
    drift_before.trackBunch(bunch)
    base_rf_gap.trackBunch(bunch)
    drift_after.trackBunch(bunch)
    delta_eKin_out_0 = bunch.getSyncParticle().kinEnergy() - Ekin_init
    exit_phase0 = (bunch.getSyncParticle().time()*cav_frequency)*180./math.pi
    exit_phase0 = phaseNearTargetPhaseDeg(exit_phase0,0.)
    
    bunch = Bunch()
    bunch_init.copyBunchTo(bunch)
    three_point_gap.trackBunch(bunch)
    delta_eKin_out_1 = bunch.getSyncParticle().kinEnergy() - Ekin_init
    exit_phase1 = (bunch.getSyncParticle().time()*cav_frequency)*180./math.pi
    exit_phase1 = phaseNearTargetPhaseDeg(exit_phase1,exit_phase0)
    
    delta_phase = exit_phase0 - exit_phase1
    
    st  = " %+8.2f  "%(cav_phase*180./math.pi)
    st += " %+8.4f  "%(delta_eKin_out_0*1000.)
    st += " %+8.4f  "%(delta_eKin_out_1*1000.)
    st += " %+8.4f  "%((delta_eKin_out_0 - delta_eKin_out_1)*1000.)
    st += "   %+8.4f   "%delta_phase 
    
    print (st)     


print ("Stop.")
sys.exit(0)