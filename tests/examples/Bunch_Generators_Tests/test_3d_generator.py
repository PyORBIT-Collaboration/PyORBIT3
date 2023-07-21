#!/usr/bin/env python

# --------------------------------------------------------
# The script will test the gerenrators from bunch_generators
# --------------------------------------------------------

import math
import random
import sys
import pytest

random.seed(100)
from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D

n = 10000
# ---------------------------------------------
# KV 3D
# ---------------------------------------------
twissX = TwissContainer(alpha=1.0, beta=2.0, emittance=3.0)
twissY = TwissContainer(alpha=2.0, beta=3.0, emittance=4.0)
twissZ = TwissContainer(alpha=3.0, beta=4.0, emittance=5.0)
dist = KVDist3D(twissX, twissY, twissZ)
twiss_analysis = TwissAnalysis(3)
for i in range(n):
    (x, xp, y, yp, z, zp) = dist.getCoordinates()
    twiss_analysis.account((x, xp, y, yp, z, zp))

print("=================================================")
print("KV 3D - done!")
print("                  alpha       beta [m/rad]    gamma     emitt[m*rad] ")
print("Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g " % twissX.getAlphaBetaGammaEmitt())
print("Generated X  %12.5g  %12.5g   %12.5g    %12.5g " % twiss_analysis.getTwiss(0))
print("......................................................................")
print("Twiss     Y  %12.5g  %12.5g   %12.5g    %12.5g " % twissY.getAlphaBetaGammaEmitt())
print("Generated Y  %12.5g  %12.5g   %12.5g    %12.5g " % twiss_analysis.getTwiss(1))
print("......................................................................")
print("Twiss     Z  %12.5g  %12.5g   %12.5g    %12.5g " % twissZ.getAlphaBetaGammaEmitt())
print("Generated Z  %12.5g  %12.5g   %12.5g    %12.5g " % twiss_analysis.getTwiss(2))
print("=================================================")

# pytest for KVDist3D
get_KV3D_Twiss_X = "%.5g %.5g %.5g %.5g" % twiss_analysis.getTwiss(0)


def test_kv3d_twissX_generator():
    assert get_KV3D_Twiss_X == "0.98779 1.9908 0.99243 2.9966"


get_KV3D_Twiss_Y = "%.5g %.5g %.5g %.5g" % twiss_analysis.getTwiss(1)


def test_kv3d_twissY_generator():
    assert get_KV3D_Twiss_Y == "2.0165 3.0098 1.6833 3.9872"


get_KV3D_Twiss_Z = "%.5g %.5g %.5g %.5g" % twiss_analysis.getTwiss(2)


def test_kv3d_twissZ_generator():
    assert get_KV3D_Twiss_Z == "3.0347 4.0329 2.5316 5.0201"


# ---------------------------------------------
# Water Bag 3D
# ---------------------------------------------
dist = WaterBagDist3D(twissX, twissY, twissZ)
twiss_analysis.init()
for i in range(n):
    (x, xp, y, yp, z, zp) = dist.getCoordinates()
    twiss_analysis.account((x, xp, y, yp, z, zp))

print("=================================================")
print("Water Bag 3D - done!")
print("                  alpha       beta [m/rad]    gamma     emitt[m*rad] ")
print("Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g " % twissX.getAlphaBetaGammaEmitt())
print("Generated X  %12.5g  %12.5g   %12.5g    %12.5g " % twiss_analysis.getTwiss(0))
print("......................................................................")
print("Twiss     Y  %12.5g  %12.5g   %12.5g    %12.5g " % twissY.getAlphaBetaGammaEmitt())
print("Generated Y  %12.5g  %12.5g   %12.5g    %12.5g " % twiss_analysis.getTwiss(1))
print("......................................................................")
print("Twiss     Z  %12.5g  %12.5g   %12.5g    %12.5g " % twissZ.getAlphaBetaGammaEmitt())
print("Generated Z  %12.5g  %12.5g   %12.5g    %12.5g " % twiss_analysis.getTwiss(2))
print("=================================================")

# Pytest for Water Bag 3D
get_WaterBag3D_Twiss_X = "%.5g %.5g %.5g %.5g" % twiss_analysis.getTwiss(0)


def test_waterbag3d_twissX_generator():
    assert get_WaterBag3D_Twiss_X == "0.98348 1.9663 1.0005 2.9886"


get_WaterBag3D_Twiss_Y = "%.5g %.5g %.5g %.5g" % twiss_analysis.getTwiss(1)


def test_waterbag3d_twissY_generator():
    assert get_WaterBag3D_Twiss_Y == "1.995 3.0085 1.6553 3.9596"


get_WaterBag3D_Twiss_Z = "%.5g %.5g %.5g %.5g" % twiss_analysis.getTwiss(2)


def test_waterbag3d_twissZ_generator():
    assert get_WaterBag3D_Twiss_Z == "3.0269 4.0174 2.5295 5.0614"


# ---------------------------------------------
# Gauss 3D
# ---------------------------------------------
dist = GaussDist3D(twissX, twissY, twissZ)
twiss_analysis.init()
for i in range(n):
    (x, xp, y, yp, z, zp) = dist.getCoordinates()
    twiss_analysis.account((x, xp, y, yp, z, zp))

print("=================================================")
print("Gauss 3D - done!")
print("                  alpha       beta [m/rad]    gamma     emitt[m*rad] ")
print("Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g " % twissX.getAlphaBetaGammaEmitt())
print("Generated X  %12.5g  %12.5g   %12.5g    %12.5g " % twiss_analysis.getTwiss(0))
print("......................................................................")
print("Twiss     Y  %12.5g  %12.5g   %12.5g    %12.5g " % twissY.getAlphaBetaGammaEmitt())
print("Generated Y  %12.5g  %12.5g   %12.5g    %12.5g " % twiss_analysis.getTwiss(1))
print("......................................................................")
print("Twiss     Z  %12.5g  %12.5g   %12.5g    %12.5g " % twissZ.getAlphaBetaGammaEmitt())
print("Generated Z  %12.5g  %12.5g   %12.5g    %12.5g " % twiss_analysis.getTwiss(2))
print("=================================================")

# pytests for GaussDist3D
get_Gauss3D_Twiss_X = "%.5g %.5g %.5g %.5g" % twiss_analysis.getTwiss(0)


def test_gauss3d_twissX_generator():
    assert get_Gauss3D_Twiss_X == "1.0149 2.0058 1.012 3.0533"


get_Gauss3D_Twiss_Y = "%.5g %.5g %.5g %.5g" % twiss_analysis.getTwiss(1)


def test_gauss3d_twissY_generator():
    assert get_Gauss3D_Twiss_Y == "2.011 3.0206 1.6699 4.0153"


get_Gauss3D_Twiss_Z = "%.5g %.5g %.5g %.5g" % twiss_analysis.getTwiss(2)


def test_gauss3d_twissZ_generator():
    assert get_Gauss3D_Twiss_Z == "3.0053 3.9921 2.5129 5.0294"
