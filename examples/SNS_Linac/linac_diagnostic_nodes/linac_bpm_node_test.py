#!/usr/bin/env python

# ------------------------------------------------------------------------
# Test of the linac BPM node. We will generate the bunch with the
# particular Twiss. The distribution function is Gaussian.
# We will test that the amplitude reported by BPM will be
# defined by the formula
# amp_bpm = const*exp(-0.5*(2*pi/360)^2*sigma2*(freq_bpm/freq_rf)^2)
# where sigma2 is a square of RMS bunch length in deg
# and in the log scale it will be
# ln(amp_bpm) = const - 0.5*(2*pi/360)^2*sigma2*(freq_bpm/freq_rf)^2
#
# At the end it will compare the phases calculated in BPM
# by different methods.
# ------------------------------------------------------------------------

import math
import sys
import os

import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

from orbit.py_linac.lattice import BaseLinacNode

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist3D
from orbit.bunch_generators import TwissAnalysis

from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit.py_linac.lattice import LinacBPM

from orbit.utils import phaseNearTargetPhaseDeg


class SNS_Linac_BunchGenerator:
    """
    Generates the pyORBIT SNS Linac Bunches using the Gauss distribution.
    Twiss parameters have the following units: x in [m], xp in [rad]
    and the X and Y emittances are un-normalized. The longitudinal emittance
    is in [GeV*m].
    """

    def __init__(self, frequency=402.5e6):
        self.bunch_frequency = frequency
        # set H- mass
        # self.bunch.mass(0.9382723 + 2*0.000511)
        self.bunch = Bunch()
        self.bunch.getSyncParticle().kinEnergy(0.0025)
        self.init_coords = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.bunch.mass(0.939294)
        self.bunch.charge(-1.0)
        self.c = 2.99792458e8  # speed of light in m/sec
        self.beam_current = 38.0  # beam current in mA , design = 38 mA
        self.rf_wave_lenght = self.c / self.bunch_frequency
        self.si_e_charge = 1.6021773e-19
        # ----------------------------------------
        self.twiss_analysis = BunchTwissAnalysis()

    def setInitialCorrdsCenter(self, x0, xp0, y0, yp0, z0, dE0):
        self.init_coords = (x0, xp0, y0, yp0, z0, dE0)

    def getInitialCorrdsCenter(self):
        return self.init_coords

    def setParticleCharge(self, charge):
        """
        Sets the particle charge H- => -1.0   and proton => +1.0
        """
        self.bunch.charge(charge)

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

    def getBunch(self, nParticles, twissX, twissY, twissZ, cut_off=-1.0):
        """
        Returns the pyORBIT bunch with particular number of particles.
        """
        (x0, xp0, y0, yp0, z0, dE0) = self.init_coords
        comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        rank = orbit_mpi.MPI_Comm_rank(comm)
        size = orbit_mpi.MPI_Comm_size(comm)
        data_type = mpi_datatype.MPI_DOUBLE
        main_rank = 0
        bunch = Bunch()
        self.bunch.copyEmptyBunchTo(bunch)
        macrosize = self.beam_current * 1.0e-3 / self.bunch_frequency
        macrosize /= math.fabs(bunch.charge()) * self.si_e_charge
        distributor = GaussDist3D(twissX, twissY, twissZ, cut_off)
        bunch.getSyncParticle().time(0.0)
        for i in range(nParticles):
            (x, xp, y, yp, z, dE) = distributor.getCoordinates()
            (x, xp, y, yp, z, dE) = orbit_mpi.MPI_Bcast((x, xp, y, yp, z, dE), data_type, main_rank, comm)
            if i % size == rank:
                bunch.addParticle(x + x0, xp + xp0, y + y0, yp + yp0, z + z0, dE + dE0)
        nParticlesGlobal = bunch.getSizeGlobal()
        bunch.macroSize(macrosize / nParticlesGlobal)
        return bunch

    def bunchAnalysis(self, bunch):
        # ---- returns (x_rms,y_rms,z_rms_deg)
        self.twiss_analysis.analyzeBunch(bunch)
        gamma = bunch.getSyncParticle().gamma()
        beta = bunch.getSyncParticle().beta()
        x_rms = math.sqrt(self.twiss_analysis.getTwiss(0)[1] * self.twiss_analysis.getTwiss(0)[3]) * 1000.0
        y_rms = math.sqrt(self.twiss_analysis.getTwiss(1)[1] * self.twiss_analysis.getTwiss(1)[3]) * 1000.0
        z_rms = math.sqrt(self.twiss_analysis.getTwiss(2)[1] * self.twiss_analysis.getTwiss(2)[3]) * 1000.0
        z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
        z_rms_deg = z_to_phase_coeff * z_rms / 1000.0
        (alphaX, betaX, emittX) = (
            self.twiss_analysis.getTwiss(0)[0],
            self.twiss_analysis.getTwiss(0)[1],
            self.twiss_analysis.getTwiss(0)[3] * 1.0e6,
        )
        (alphaY, betaY, emittY) = (
            self.twiss_analysis.getTwiss(1)[0],
            self.twiss_analysis.getTwiss(1)[1],
            self.twiss_analysis.getTwiss(1)[3] * 1.0e6,
        )
        (alphaZ, betaZ, emittZ) = (
            self.twiss_analysis.getTwiss(2)[0],
            self.twiss_analysis.getTwiss(2)[1],
            self.twiss_analysis.getTwiss(2)[3] * 1.0e6,
        )
        norm_emittX = emittX * gamma * beta
        norm_emittY = emittY * gamma * beta
        # gammaZ = (1.0+alphaZ**2)/betaZ
        # dE_sigma = math.sqrt(gammaZ*emittZ)
        # print "debug dE_sigma [MeV]=",dE_sigma
        # ---- phi_de_emittZ will be in [pi*deg*MeV]
        phi_de_emittZ = z_to_phase_coeff * emittZ
        return (x_rms, y_rms, z_rms_deg)

    def shiftPhase(self, bunch, delta_phi_deg):
        dz = -delta_phi_deg / self.getZtoPhaseCoeff(bunch)
        nParts = bunch.getSize()
        for ind in range(nParts):
            z = bunch.z(ind)
            bunch.z(ind, z + dz)


# -------------------------------------------
# START of Script
# -------------------------------------------

comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)
size = orbit_mpi.MPI_Comm_size(comm)
data_type = mpi_datatype.MPI_DOUBLE
main_rank = 0

(alphaX, betaX, emittX) = (3.25, 1.95, 3.37 * 1.0e-6)
(alphaY, betaY, emittY) = (1.18, 0.42, 3.09 * 1.0e-6)
(alphaZ, betaZ, emittZ) = (-0.82, 232.2, 0.0165 * 1.0e-6)

print(" ========= PyORBIT Twiss ===========")
print(" aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f " % (alphaX, betaX, emittX * 1.0e6))
print(" aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f " % (alphaY, betaY, emittY * 1.0e6))
print(" aplha beta emitt[mm*MeV] Z= %6.4f %6.4f %6.4f " % (alphaZ, betaZ, emittZ * 1.0e6))

twissX = TwissContainer(alphaX, betaX, emittX)
twissY = TwissContainer(alphaY, betaY, emittY)
twissZ = TwissContainer(alphaZ, betaZ, emittZ)

print("Start Bunch Generation.")
bunch_gen = SNS_Linac_BunchGenerator()

nParts = 100000

bunch = bunch_gen.getBunch(nParts, twissX, twissY, twissZ)
(x_rms, y_rms, z_rms_deg) = bunch_gen.bunchAnalysis(bunch)
print("init bunch sizes = %4.2f %4.2f %4.2f " % (x_rms, y_rms, z_rms_deg))

z_rms_deg_ini = z_rms_deg
max_z_rms_deg = 100.0
z_rms_deg_step = 5.0
z_rms_deg_new = 0.0

bpm_node = LinacBPM(frequency=805.0e6, name="BPM")

# -------------------------------------------
#  Check the BPM amplitude calculations
# -------------------------------------------
print("========================================================================")
print(" size(deg)   size^2  bpm_amp   ln(bpm_amp)  ln(bpm_amp)_theory   diff ")

while z_rms_deg_new < max_z_rms_deg:
    z_rms_deg_new += z_rms_deg_step
    emittZ_new = emittZ * (z_rms_deg_new / z_rms_deg_ini) ** 2
    twissZ = TwissContainer(alphaZ, betaZ, emittZ_new)
    bunch = bunch_gen.getBunch(nParts, twissX, twissY, twissZ)
    (x_rms, y_rms, z_rms_deg) = bunch_gen.bunchAnalysis(bunch)
    bpm_node.analyzeBunch(bunch)
    amp = bpm_node.getAmplitude()
    ln_amp = math.log(amp)
    y_theory = -2 * (2 * math.pi / 360.0) ** 2 * z_rms_deg**2
    diff = -(ln_amp - y_theory)
    print(" %6.4f  %8.4f  %12.5g    %12.5g  %12.5g   %12.5g " % (z_rms_deg, z_rms_deg**2, amp, ln_amp, y_theory, diff))

# -------------------------------------------
#  Check the BPM phase calculations
# -------------------------------------------
twissZ = TwissContainer(alphaZ, betaZ, emittZ)
bunch = bunch_gen.getBunch(nParts, twissX, twissY, twissZ)

# ---- the real BPM phase shift will be 20 deg because the bunch ferquency (RF)
# ---- is two times less than BPM frequency
time_shift = (1.0 / 402.5e6) * (10.0 / 360.0)
bunch.getSyncParticle().time(time_shift)

phase_avg_old = 0.0
phase_peak_old = 0.0
phase_fourier_old = 0.0
phase_synch_old = 0.0
phase_rms_old = 0.0

phase_step = 5.0
phase = 0.0

print("========================================================================")
print(" phase amp   phase_rms   phase_synch   phase_avg   phase_peak   phase_fourier ")
for phase_ind in range(int(360 / phase_step)):
    phase += phase_step
    bunch_gen.shiftPhase(bunch, phase_step)
    bpm_node.analyzeBunch(bunch)
    phase_avg = phaseNearTargetPhaseDeg(bpm_node.getAvgPhase(), phase_avg_old)
    phase_peak = phaseNearTargetPhaseDeg(bpm_node.getPeakPhase(), phase_peak_old)
    phase_fourier = phaseNearTargetPhaseDeg(bpm_node.getFourierPhase(), phase_fourier_old)
    phase_synch = phaseNearTargetPhaseDeg(bpm_node.getSynchPhase(), phase_synch_old)
    phase_rms = phaseNearTargetPhaseDeg(bpm_node.getPhaseRMS(), phase_rms_old)
    amp = bpm_node.getAmplitude()
    phase_avg_old = phase_avg
    phase_peak_old = phase_peak
    phase_fourier_old = phase_fourier
    phase_synch_old = phase_synch
    phase_rms_old = phase_rms
    st = " %6.2f " % phase
    st += "  %6.5f  " % amp
    st += "   %6.1f    %6.1f  %6.1f  %6.1f  %6.1f  " % (phase_rms, phase_synch, phase_avg, phase_peak, phase_fourier)
    print(st)


sys.exit(0)

# ---------------------------------------------------------
#   Results (remember about randomness in the bunch generation)
# ---------------------------------------------------------
"""[shi@shilnx linac_diagnostics_nodes]$ pyORBIT linac_bpm_node_test.py
 ========= PyORBIT Twiss ===========
 aplha beta emitt[mm*mrad] X= 3.2500 1.9500 3.3700
 aplha beta emitt[mm*mrad] Y= 1.1800 0.4200 3.0900
 aplha beta emitt[mm*MeV] Z= -0.8200 232.2000 0.0165
Start Bunch Generation.
init bunch sizes = 2.56 1.13 12.98
========================================================================
 size(deg)   size^2  bpm_amp   ln(bpm_amp)  ln(bpm_amp)_theory   diff
 5.0134   25.1339       0.31347         -1.1601     -0.015312         1.1447
 10.0184  100.3683       0.29942         -1.2059     -0.061148         1.1448
 15.0564  226.6951       0.27724         -1.2829      -0.13811         1.1448
 20.0492  401.9709       0.24913         -1.3898      -0.24489         1.1449
 25.0748  628.7473       0.21705         -1.5276      -0.38305         1.1446
 30.1575  909.4737       0.18302         -1.6982      -0.55408         1.1441
 34.9025  1218.1849       0.15153          -1.887      -0.74216         1.1448
 40.2496  1620.0284       0.11892         -2.1293      -0.98698         1.1423
 44.9608  2021.4734      0.092611         -2.3793       -1.2316         1.1478
 50.0234  2502.3437      0.069563         -2.6655       -1.5245          1.141
 54.9914  3024.0567       0.04999         -2.9959       -1.8424         1.1536
 60.1796  3621.5787      0.035455         -3.3395       -2.2064         1.1331
 65.3219  4266.9513      0.023545         -3.7488       -2.5996         1.1493
 69.9689  4895.6426      0.016286         -4.1175       -2.9826         1.1349
 75.2341  5660.1759      0.010361         -4.5697       -3.4484         1.1213
 80.3607  6457.8438     0.0062813         -5.0702       -3.9343         1.1358
 85.0716  7237.1832     0.0041836         -5.4766       -4.4091         1.0674
 89.8972  8081.5040      0.002281         -6.0831       -4.9235         1.1596
 95.1672  9056.7999     0.0026778         -5.9228       -5.5177        0.40505
 100.1021  10020.4334    0.00098177         -6.9262       -6.1048        0.82136
========================================================================
 phase amp   phase_rms   phase_synch   phase_avg   phase_peak   phase_fourier
   5.00   0.28709       26.0      20.0    30.1    26.5    30.1
  10.00   0.28709       26.0      20.0    40.1    36.5    40.1
  15.00   0.28709       26.0      20.0    50.1    46.5    50.1
  20.00   0.28709       26.0      20.0    60.1    56.5    60.1
  25.00   0.28709       26.0      20.0    70.1    66.5    70.1
  30.00   0.28709       26.0      20.0    80.1    76.5    80.1
  35.00   0.28709       26.0      20.0    90.1    86.5    90.1
  40.00   0.28709       26.0      20.0   100.1    96.5   100.1
  45.00   0.28709       26.0      20.0   110.1   106.5   110.1
  50.00   0.28709       26.0      20.0   120.1   116.5   120.1
  55.00   0.28709       26.0      20.0   130.1   126.5   130.1
  60.00   0.28709       26.0      20.0   140.1   136.5   140.1
  65.00   0.28709       26.0      20.0   150.1   146.5   150.1
  70.00   0.28709       26.0      20.0   160.1   156.5   160.1
  75.00   0.28709       26.0      20.0   170.1   166.5   170.1
  80.00   0.28709       26.0      20.0   180.1   176.5   180.1
  85.00   0.28709       26.0      20.0   190.1   186.5   190.1
  90.00   0.28709       26.0      20.0   200.1   196.5   200.1
  95.00   0.28709       26.0      20.0   210.1   206.5   210.1
 100.00   0.28709       26.0      20.0   220.1   216.5   220.1
 105.00   0.28709       26.0      20.0   230.1   226.5   230.1
 110.00   0.28709       26.0      20.0   240.1   236.5   240.1
 115.00   0.28709       26.0      20.0   250.1   246.5   250.1
 120.00   0.28709       26.0      20.0   260.1   256.5   260.1
 125.00   0.28709       26.0      20.0   270.1   266.5   270.1
 130.00   0.28709       26.0      20.0   280.1   276.5   280.1
 135.00   0.28709       26.0      20.0   290.1   286.5   290.1
 140.00   0.28709       26.0      20.0   300.1   296.5   300.1
 145.00   0.28709       26.0      20.0   310.1   306.5   310.1
 150.00   0.28709       26.0      20.0   320.1   316.5   320.1
 155.00   0.28709       26.0      20.0   330.1   326.5   330.1
 160.00   0.28709       26.0      20.0   340.1   336.5   340.1
 165.00   0.28709       26.0      20.0   350.1   346.5   350.1
 170.00   0.28709       26.0      20.0   360.1   356.5   360.1
 175.00   0.28709       26.0      20.0   370.1   366.5   370.1
 180.00   0.28709       26.0      20.0   380.1   376.5   380.1
 185.00   0.28709       26.0      20.0   390.1   386.5   390.1
 190.00   0.28709       26.0      20.0   400.1   396.5   400.1
 195.00   0.28709       26.0      20.0   410.1   406.5   410.1
 200.00   0.28709       26.0      20.0   420.1   416.5   420.1
 205.00   0.28709       26.0      20.0   430.1   426.5   430.1
 210.00   0.28709       26.0      20.0   440.1   436.5   440.1
 215.00   0.28709       26.0      20.0   450.1   446.5   450.1
 220.00   0.28709       26.0      20.0   460.1   456.5   460.1
 225.00   0.28709       26.0      20.0   470.1   466.5   470.1
 230.00   0.28709       26.0      20.0   480.1   476.5   480.1
 235.00   0.28709       26.0      20.0   490.1   486.5   490.1
 240.00   0.28709       26.0      20.0   500.1   496.5   500.1
 245.00   0.28709       26.0      20.0   510.1   506.5   510.1
 250.00   0.28709       26.0      20.0   520.1   516.5   520.1
 255.00   0.28709       26.0      20.0   530.1   526.5   530.1
 260.00   0.28709       26.0      20.0   540.1   536.5   540.1
 265.00   0.28709       26.0      20.0   550.1   546.5   550.1
 270.00   0.28709       26.0      20.0   560.1   556.5   560.1
 275.00   0.28709       26.0      20.0   570.1   566.5   570.1
 280.00   0.28709       26.0      20.0   580.1   576.5   580.1
 285.00   0.28709       26.0      20.0   590.1   586.5   590.1
 290.00   0.28709       26.0      20.0   600.1   596.5   600.1
 295.00   0.28709       26.0      20.0   610.1   606.5   610.1
 300.00   0.28709       26.0      20.0   620.1   616.5   620.1
 305.00   0.28709       26.0      20.0   630.1   626.5   630.1
 310.00   0.28709       26.0      20.0   640.1   636.5   640.1
 315.00   0.28709       26.0      20.0   650.1   646.5   650.1
 320.00   0.28709       26.0      20.0   660.1   656.5   660.1
 325.00   0.28709       26.0      20.0   670.1   666.5   670.1
 330.00   0.28709       26.0      20.0   680.1   676.5   680.1
 335.00   0.28709       26.0      20.0   690.1   686.5   690.1
 340.00   0.28709       26.0      20.0   700.1   696.5   700.1
 345.00   0.28709       26.0      20.0   710.1   706.5   710.1
 350.00   0.28709       26.0      20.0   720.1   716.5   720.1
 355.00   0.28709       26.0      20.0   730.1   726.5   730.1
 360.00   0.28709       26.0      20.0   740.1   736.5   740.1
"""
