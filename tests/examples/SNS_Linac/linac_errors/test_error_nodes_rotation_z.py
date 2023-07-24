import pytest
import math
import sys
import os
import random
import orbit.core

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer

from orbit.py_linac.lattice import LinacAccLattice
from orbit.py_linac.lattice import Sequence
from orbit.py_linac.lattice import Drift
from orbit.py_linac.lattice import Quad
from orbit.py_linac.lattice import Bend

from orbit.py_linac.errors import ErrorCntrlStraightRotationX
from orbit.py_linac.errors import ErrorCntrlStraightRotationY
from orbit.py_linac.errors import ErrorCntrlStraightRotationZ

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist3D
from orbit.bunch_generators import TwissAnalysis

from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit_utils import Matrix, PhaseVector

random.seed(100)


def printMatrix(m, name=""):
    """
    Prints matrix
    """
    print("----matrix--- size=", m.size(), " name=", name)
    for i in range(m.size()[0]):
        for j in range(m.size()[1]):
            st = "m(%1d,%1d) = %+12.5g " % (i, j, m.get(i, j))
            print(st, end=" ")
        print("")


def matrixToString(m, name=""):
    """
    returns matrix as a string
    """
    matrix_string = ""
    # matrix_string = "----matrix--- size= {}  name= {}\n".format(m.size(), name)
    for i in range(m.size()[0]):
        for j in range(m.size()[1]):
            st = "m(%1d,%1d) = %+12.5g " % (i, j, m.get(i, j))
            matrix_string += st + " "
        matrix_string += "\n"
    return matrix_string


def setM(arr):
    """
    Sets up matrix from 2D array
    """
    n = len(arr)
    m = Matrix(n, n)
    for i in range(n):
        for j in range(n):
            m.set(i, j, arr[i][j])
    return m


def getQuadMatrix(bunch, gradient, length):
    """
    Returns transport matrix for quad
    """
    momentum = bunch.getSyncParticle().momentum()
    charge = bunch.charge()
    kqc = charge * gradient / (3.335640952 * momentum)
    sqrt_kq = math.sqrt(abs(kqc))
    kqlength = sqrt_kq * length
    cn = math.cos(kqlength)
    sn = math.sin(kqlength)
    ch = math.cosh(kqlength)
    sh = math.sinh(kqlength)
    arr = [[], [], [], []]
    if kqc > 0.0:
        arr[0] = [cn, sn / sqrt_kq, 0.0, 0.0]
        arr[1] = [-sn * sqrt_kq, cn, 0.0, 0.0]
        arr[2] = [0.0, 0.0, ch, sh / sqrt_kq]
        arr[3] = [0.0, 0.0, sh * sqrt_kq, ch]
    else:
        arr[0] = [ch, sh / sqrt_kq, 0.0, 0.0]
        arr[1] = [sh * sqrt_kq, ch, 0.0, 0.0]
        arr[2] = [0.0, 0.0, cn, sn / sqrt_kq]
        arr[3] = [0.0, 0.0, -sn * sqrt_kq, cn]
    return setM(arr)


def getRotationMatrixZ(angle):
    """
    Returns rotation matrix for Z-axis rotation
    """
    cs = math.cos(angle)
    sn = math.sin(angle)
    arr = [[], [], [], []]
    arr[0] = [cs, 0.0, sn, 0.0]
    arr[1] = [0.0, cs, 0.0, sn]
    arr[2] = [-sn, 0.0, cs, 0.0]
    arr[3] = [0.0, -sn, 0.0, cs]
    return setM(arr)


def trackBunch4D(bunch, matrix):
    """
    Tracks bunch particles through the 4D matrix for x,xp,y,yp.
    Returns new bunch with new coordinates.
    """
    m = matrix
    bunch_out = Bunch()
    bunch.copyEmptyBunchTo(bunch_out)
    vct = PhaseVector(4)
    for ind in range(bunch.getSize()):
        vct.set(0, bunch.x(ind))
        vct.set(1, bunch.xp(ind))
        vct.set(2, bunch.y(ind))
        vct.set(3, bunch.yp(ind))
        # -----------------------------
        vout = matrix.mult(vct)
        # -----------------------------
        z = bunch.z(ind)
        dE = bunch.dE(ind)
        bunch_out.addParticle(vout.get(0), vout.get(1), vout.get(2), vout.get(3), z, dE)
    return bunch_out


def getBunch4D_CorrelationMatrix(twiss_analysis, bunch):
    """
    Returns the 4D correlation matrix
    """
    twiss_analysis.analyzeBunch(bunch)
    arr = []
    sigma_arr = []
    for row_ind in range(4):
        arr.append([])
        for col_ind in range(4):
            corr = twiss_analysis.getCorrelation(row_ind, col_ind)
            arr[row_ind].append(corr)
            if row_ind == col_ind:
                sigma_arr.append(math.sqrt(corr))
    for row_ind in range(4):
        for col_ind in range(4):
            arr[row_ind][col_ind] /= sigma_arr[row_ind] * sigma_arr[col_ind]
    corr_matrix = setM(arr)
    return corr_matrix


def bunchCentering(bunch):
    """
    Bunch center after generating can have small deviation from the (0,0,0,0,0,0)
    This function will force centering the bunch.
    """
    twiss_analysis = BunchTwissAnalysis()
    twiss_analysis.analyzeBunch(bunch)
    # -----------------------------------------------
    # let's center the beam
    (x_avg, y_avg) = (twiss_analysis.getAverage(0), twiss_analysis.getAverage(2))
    (xp_avg, yp_avg) = (twiss_analysis.getAverage(1), twiss_analysis.getAverage(3))
    (z_avg, dE_avg) = (twiss_analysis.getAverage(4), twiss_analysis.getAverage(5))
    for part_id in range(bunch.getSize()):
        bunch.x(part_id, bunch.x(part_id) - x_avg)
        bunch.y(part_id, bunch.y(part_id) - y_avg)
        bunch.xp(part_id, bunch.xp(part_id) - xp_avg)
        bunch.yp(part_id, bunch.yp(part_id) - yp_avg)
        bunch.z(part_id, bunch.z(part_id) - z_avg)
        bunch.dE(part_id, bunch.dE(part_id) - dE_avg)
    # -----------------------------------------------
    return (x_avg, y_avg, dE_avg)


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
        self.init_coords = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.bunch.mass(0.939294)
        self.bunch.charge(-1.0)
        self.bunch.getSyncParticle().kinEnergy(0.0025)
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
        bunch = Bunch()
        self.bunch.copyEmptyBunchTo(bunch)
        macrosize = self.beam_current * 1.0e-3 / self.bunch_frequency
        macrosize /= math.fabs(bunch.charge()) * self.si_e_charge
        distributor = GaussDist3D(twissX, twissY, twissZ, cut_off)
        bunch.getSyncParticle().time(0.0)
        for i in range(nParticles):
            (x, xp, y, yp, z, dE) = distributor.getCoordinates()
            bunch.addParticle(x + x0, xp + xp0, y + y0, yp + yp0, z + z0, dE + dE0)
        nParticlesGlobal = bunch.getSizeGlobal()
        bunch.macroSize(macrosize / nParticlesGlobal)
        return bunch


# ---------------------------------------------------------------
# ----    START of the script
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# ---- Let's make a linac lattice
# ---------------------------------------------------------------

drift_1 = Drift("drift_1")
drift_2 = Drift("drift_2")

drift_1.setLength(0.5)
drift_2.setLength(0.5)

drift_1.setnParts(10)
drift_2.setnParts(10)

quad_1 = Quad("quad_1")
quad_1.setLength(0.05)
gradient = -5.0
quad_1.setParam("dB/dr", gradient)
quad_1.setnParts(6)

accSeq = Sequence("Test_3_Drifts_AccSeq")
accSeq.addNode(drift_1)
accSeq.addNode(quad_1)
accSeq.addNode(drift_2)

accLattice = LinacAccLattice("Error_3_Drifts_Lattice")

for node in accSeq.getNodes():
    accLattice.addNode(node)
accLattice.initialize()

# ---- This will force the usage of linac type quad tracker
# ---- This tracker is defined in /scr/linac/tracking/linac_tracking.cc
# ---- If it is not True, the usual Teapot like quad tracking functions
# ---- will be used from /scr/teapot/teapotbase.cc
accLattice.setLinacTracker(True)

print("==============Lattice======================")
nodes = accLattice.getNodes()
for node in nodes:
    print("node =", node.getName(), " length=", node.getLength())
print("quad gradient [T/m] = ", gradient)

rotation_angle = math.pi / 4
print("quad rotation angle [grad] = %6.3f " % (rotation_angle * 180.0 / math.pi))

print("===========================================")
# ---------------------------------------------------------------------------
# Here the example how to assign the error functions to a particular node
# in the lattice.
# ---------------------------------------------------------------------------
node_with_error = accLattice.getNodeForName("quad_1")

errorCntrl_z = ErrorCntrlStraightRotationZ(node_with_error.getName() + "_errorNodeCntrl")
errorCntrl_z.setOneNodeParent(node_with_error)
errorCntrl_z.setRotationAngle(rotation_angle)  # rotation angle in radians

# -------------------------------------------------------
#    Now let's generate bunch - nothing new
# -------------------------------------------------------

bunch_generator = SNS_Linac_BunchGenerator()

peak_current = 0.0  # mA
bunch_generator.setBeamCurrent(peak_current)

# ------ PyORBIT emittances
(alphaX, betaX, emittX) = (-1.0, 4.0, 1.0 * 1.0e-6)
(alphaY, betaY, emittY) = (+1.0, 2.0, 1.0 * 1.0e-6)
(alphaZ, betaZ, emittZ) = (-0.02, 100.0, 0.016 * 1.0e-6)

print(" ========= PyORBIT Twiss ===========")
print(" aplha beta emitt[mm*mrad] X= %+6.4f %6.4f %6.4f " % (alphaX, betaX, emittX * 1.0e6))
print(" aplha beta emitt[mm*mrad] Y= %+6.4f %6.4f %6.4f " % (alphaY, betaY, emittY * 1.0e6))
print(" aplha beta emitt[m*GeV]   Z= %+6.4f %6.2f %6.4f " % (alphaZ, betaZ, emittZ * 1.0e6))
print("=================================================")

twissX = TwissContainer(alphaX, betaX, emittX)
twissY = TwissContainer(alphaY, betaY, emittY)
twissZ = TwissContainer(alphaZ, betaZ, emittZ)

nParticles = 50000

bunch_init = bunch_generator.getBunch(nParticles, twissX, twissY, twissZ)

# ---- center the bunch particle distribution
(x_avg, y_avg, dE_avg) = bunchCentering(bunch_init)

bunch = Bunch()
bunch_init.copyBunchTo(bunch)

twiss_analysis = BunchTwissAnalysis()

corr_matrix = getBunch4D_CorrelationMatrix(twiss_analysis, bunch)
printMatrix(corr_matrix, " Initial Correlation Matrix")

accLattice.trackDesignBunch(bunch)
accLattice.trackBunch(bunch)

print("=================================================")
corr_matrix = getBunch4D_CorrelationMatrix(twiss_analysis, bunch)
printMatrix(corr_matrix, " Final Correlation Matrix by TEAPOT tracking")
print("=================================================")
corr_matrix_teapot_track = matrixToString(corr_matrix, " Final Correlation Matrix by TEAPOT tracking")
L = drift_1.getLength()
arr = [[], [], [], []]
arr[0] = [1.0, L, 0.0, 0.0]
arr[1] = [0.0, 1.0, 0.0, 0.0]
arr[2] = [0.0, 0.0, 1.0, L]
arr[3] = [0.0, 0.0, 0.0, 1.0]
drift_1_matrix = setM(arr)
# printMatrix(drift_1_matrix,drift_1.getName())

L = quad_1.getLength()
G = quad_1.getParam("dB/dr")
quad_1_matrix = getQuadMatrix(bunch, G, L)
# printMatrix(quad_1_matrix,quad_1.getName())

L = drift_2.getLength()
arr = [[], [], [], []]
arr[0] = [1.0, L, 0.0, 0.0]
arr[1] = [0.0, 1.0, 0.0, 0.0]
arr[2] = [0.0, 0.0, 1.0, L]
arr[3] = [0.0, 0.0, 0.0, 1.0]
drift_2_matrix = setM(arr)
# printMatrix(drift_2_matrix,drift_2.getName())

rotationMtrx1 = getRotationMatrixZ(rotation_angle)
rotationMtrx2 = getRotationMatrixZ(-rotation_angle)

# ---- let's calculate the total trasport matrix for the lattice
trMatrixTotal = drift_2_matrix.mult(rotationMtrx2)
trMatrixTotal = trMatrixTotal.mult(quad_1_matrix)
trMatrixTotal = trMatrixTotal.mult(rotationMtrx1)
trMatrixTotal = trMatrixTotal.mult(drift_1_matrix)

print("=================================================")

# printMatrix(trMatrixTotal," Lattice matrix")

# ---- track bunch through the total transport matrix
bunch = trackBunch4D(bunch_init, trMatrixTotal)
corr_matrix = getBunch4D_CorrelationMatrix(twiss_analysis, bunch)
printMatrix(corr_matrix, " Final Correlation Matrix by using linear transport model")

corr_matrix_linear_track = matrixToString(corr_matrix, " Final Correlation Matrix by using linear transport model")


def test_teapotMatrix():
    expected = """m(0,0) =           +1  m(0,1) =     +0.62698  m(0,2) =     -0.89476  m(0,3) =     -0.85977  
m(1,0) =     +0.62698  m(1,1) =           +1  m(1,2) =     -0.78602  m(1,3) =     -0.24323  
m(2,0) =     -0.89476  m(2,1) =     -0.78602  m(2,2) =           +1  m(2,3) =     +0.74658  
m(3,0) =     -0.85977  m(3,1) =     -0.24323  m(3,2) =     +0.74658  m(3,3) =           +1  \n"""

    assert corr_matrix_teapot_track == expected


def test_corrMatrix():
    expected = """m(0,0) =           +1  m(0,1) =     +0.62698  m(0,2) =     -0.89476  m(0,3) =     -0.85977  
m(1,0) =     +0.62698  m(1,1) =           +1  m(1,2) =     -0.78602  m(1,3) =     -0.24323  
m(2,0) =     -0.89476  m(2,1) =     -0.78602  m(2,2) =           +1  m(2,3) =     +0.74658  
m(3,0) =     -0.85977  m(3,1) =     -0.24323  m(3,2) =     +0.74658  m(3,3) =           +1  \n"""

    assert corr_matrix_linear_track == expected
