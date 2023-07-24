# ------------------------------------------------------
# This is an example of Field Sources of electric and magnetic
# fields for 3D tracking using RK4_Tracker
# -------------------------------------------------------

import sys
import math
import time
import orbit.core
import pytest

# from orbit_utils import field_sources
from field_sources import QuadFieldSource
from field_sources import DipoleFieldSource
from field_sources import MagnetFieldSourceGrid3D

from orbit_utils import Matrix, PhaseVector
from spacecharge import Grid3D

transfCoordsMatrix = Matrix(4, 4)
transfCoordsMatrix.unit()

# ----------------------------------------------
#   Quad set get parameters test
# ----------------------------------------------

quad_field_source = QuadFieldSource()
quad_field_source.length(1.5)
print("quad length =", quad_field_source.length())
quad_field_source.gradient(3.0)
print("quad gradient =", quad_field_source.gradient())

quad_field_source.transormfMatrix(transfCoordsMatrix)
transfCoordsMatrix = quad_field_source.transormfMatrix()
transfCoordsMatrix.unit()

# ----------------------------------------------
#   Dipole set get parameters test
# ----------------------------------------------

dipole_field_source = DipoleFieldSource()
dipole_field_source.fieldsXYZ(0.1, 0.2, 0.3)
print("dipole fieldsXYZ = ", dipole_field_source.fieldsXYZ())
dipole_field_source.sizesXYZ(0.5, 0.6, 0.7)
print("dipole sizesXYZ = ", dipole_field_source.sizesXYZ())

dipole_field_source.transormfMatrix(transfCoordsMatrix)
transfCoordsMatrix = dipole_field_source.transormfMatrix()
transfCoordsMatrix.unit()

# ----------------------------------------------
#   Grid3D field source set get parameters test
# ----------------------------------------------

bx_grid3d = Grid3D(10, 10, 10)
by_grid3d = Grid3D(10, 10, 10)
bz_grid3d = Grid3D(10, 10, 10)

magnet_field3D = MagnetFieldSourceGrid3D(bx_grid3d, by_grid3d, bz_grid3d)

magnet_field3D.transormfMatrix(transfCoordsMatrix)
transfCoordsMatrix = magnet_field3D.transormfMatrix()
transfCoordsMatrix.unit()

print("Stop")


def test_quad_field_source_length():
    assert quad_field_source.length() == 1.5


def test_quad_field_source_gradient():
    assert quad_field_source.gradient() == 3.0


def test_dipole_field_sources_fieldsXYZ():
    assert dipole_field_source.fieldsXYZ() == (0.1, 0.2, 0.3)


def test_dipole_field_sources_sizesXYZ():
    assert dipole_field_source.sizesXYZ() == (0.5, 0.6, 0.7)
