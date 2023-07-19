import orbit.core
import pytest

from orbit.core import _orbit

from orbit.py_linac.lattice import Quad
from orbit.teapot import QuadTEAPOT

from bunch import Bunch


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


myfile = "out.txt"
myfile2 = "out2.txt"

bunch = Bunch()
syncPart = bunch.getSyncParticle()
syncPart.kinEnergy(1.0)
(x, xp, y, yp, z, dE) = (0.1, 0.1, 0.1, 0.1, 0.1, 0.001)
bunch.addParticle(x, xp, y, yp, z, dE)

# bunch.dumpBunch("out.txt")

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
bunch_linac.dumpBunch(myfile)

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

bunch_teapot.dumpBunch(myfile2)


str = read_lines(myfile)
print(str)
str2 = read_lines(myfile2)
print(str2)


# pytest things
@pytest.fixture
def output():
    return str


def test_multipole_example_output_linacquad(output):
    assert output == "0.14635674 0.10095515 0.1879077 0.27735323 0.089551404 0.001"


@pytest.fixture
def output2():
    return str2


def test_multipole_example_output_teapotquad(output2):
    assert output2 == "0.14635469 0.10098071 0.18790445 0.27736948 0.089552597 0.001"
