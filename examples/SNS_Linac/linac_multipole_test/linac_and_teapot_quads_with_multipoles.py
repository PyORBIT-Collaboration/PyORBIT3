#! /usr/bin/env python

"""
This script will track the bunch with one macroparticle
through the linac and TEAPOT quad to compare the results for effects
from multipole components.

The linac integration scheme is a 3-point TEAPOT scheme applied
several times.

Increase of the number of parts for each quad improve results.

"""

import sys
import math
import random
import time

from orbit.py_linac.lattice import Quad

from orbit.teapot import QuadTEAPOT

from bunch import Bunch

# ------------------------------------------
# Proton bunch with 1 GeV energy
# ------------------------------------------
bunch = Bunch()
syncPart = bunch.getSyncParticle()
syncPart.kinEnergy(1.0)
(x, xp, y, yp, z, dE) = (0.1, 0.1, 0.1, 0.1, 0.1, 0.001)
bunch.addParticle(x, xp, y, yp, z, dE)

# ------------------------------------------
# Two quads: linac and TEAPOT types
# ------------------------------------------
linac_quad = Quad("Linac_Quad")
teapot_quad = QuadTEAPOT("TEAPOT_Quad")

# ---- make copies of bunch for linac and TEAPOT quads
bunch_linac = Bunch()
bunch_teapot = Bunch()
bunch.copyBunchTo(bunch_linac)
bunch.copyBunchTo(bunch_teapot)

# --------------------------------------------
# ---- LINAC Quad    start
# --------------------------------------------
n_quad_parts = 10
G = 5.0  # T/m
paramsDict = {}
paramsDict["bunch"] = bunch_linac
linac_quad.addParam("dB/dr", G)

# ---- multipole components
# ---- poleArr = [2,3] - sextupoles and octupoles
poleArr = [2, 3]
klArr = [5.0, 50.0]
skewArr = [0, 0, 0]
linac_quad.setParam("poles", poleArr)
linac_quad.setParam("kls", klArr)
linac_quad.setParam("skews", skewArr)

linac_quad.setnParts(n_quad_parts)
linac_quad.setLength(0.5)
linac_quad.initialize()
n_parts = linac_quad.getnParts()
for ind in range(n_parts):
    linac_quad.setActivePartIndex(ind)
    linac_quad.track(paramsDict)
print("================ Linac Quad Tracking =============")
bunch_linac.dumpBunch()
print("==================================================")


# --------------------------------------------
# ---- TEAPOT Quad    start
# --------------------------------------------
n_quad_parts = 30
momentum = bunch_teapot.getSyncParticle().momentum()
kq = linac_quad.getParam("dB/dr") / (3.335640952 * momentum)
paramsDict["bunch"] = bunch_teapot
teapot_quad.addParam("kq", kq)

# ---- multipole components
# ---- poleArr = [2,3] - sextupoles and octupoles
poleArr = [2, 3]
klArr = [5.0, 50.0]
skewArr = [0, 0, 0]
teapot_quad.setParam("poles", poleArr)
teapot_quad.setParam("kls", klArr)
teapot_quad.setParam("skews", skewArr)

teapot_quad.setnParts(n_quad_parts)
teapot_quad.setLength(0.5)
teapot_quad.initialize()
n_parts = teapot_quad.getnParts()
for ind in range(n_parts):
    teapot_quad.setActivePartIndex(ind)
    teapot_quad.track(paramsDict)
print("================ TEAPOT Quad Tracking ============")
bunch_teapot.dumpBunch()
print("==================================================")
