# ------------------------------------------------------
# This is an example of Function and SplineCH
# They are containers for (x,y) (or (x,y,err) ) points
# Function provides a linear interpolation, and SplineCH
# uses 3-rd order polynomials. SplineCH can be used for
# derivatives calculations.
# -------------------------------------------------------

import sys
import math
import orbit.core
import pytest

from orbit_utils import Function
from orbit_utils import SplineCH

f = Function()


def convert_to_wanted(number):
    return f"{number:.12g}"


def FF(x):
    return math.sin(x)


def FFP(x):
    return math.cos(x)


n = 40
step = 2 * math.pi / n
for i in range(n):
    x = step * i + 0.1 * ((1.0 * i) / n) ** 2
    y = FF(x)
    f.add(x, y)

f.dump()

spline = SplineCH()
spline.compile(f)

spline.dump()

print("================")
n = 100
step = 0.8 * (f.getMaxX() - f.getMinX()) / (n - 1)
y_dev_max = 0.0
yp_dev_max = 0.0
for j in range(n):
    x = f.getMinX() + j * step
    y = FF(x)
    yp = FFP(x)
    ys = math.fabs(spline.getY(x) - y)
    yps = math.fabs(spline.getYP(x) - yp)
    if y_dev_max < ys:
        # print "debug x=",x," dev y=",ys," y=",y
        y_dev_max = ys
    if yp_dev_max < yps:
        # print "debug x=",x," dev yp=",yps," yp=",yp
        yp_dev_max = yps

print(" deviation y max =", y_dev_max)
print(" deviation yp max =", yp_dev_max)


@pytest.fixture
def get_y_dev_max():
    wantedd = convert_to_wanted(y_dev_max)
    return wantedd


def test_deviation_y_max(get_y_dev_max):
    assert get_y_dev_max == "6.34709129114e-05"


@pytest.fixture
def get_yp_dev_max():
    wanted = convert_to_wanted(yp_dev_max)
    return wanted


def test_deviation_yp_max(get_yp_dev_max):
    assert get_yp_dev_max == "0.00415191930103"
