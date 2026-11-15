#!/bin/bash

set -e -x

python fodo.py --sc 0
python fodo.py --sc 1
python fodo.py --sc 1 --offset-x 0.001
python fodo.py --sc 1 --tilt 45.0
python fodo_match.py
python fodo_optimize.py
python fodo_speed.py --sc 0
python fodo_speed.py --sc 1

python uniform_ellipsoid.py
python uniform_ellipsoid_pic.py --sc 0
python uniform_ellipsoid_pic.py --sc 1
python uniform_ellipsoid_pic.py --sc 1 --rms-y 0.002 --rot-z 45.0
python uniform_ellipsoid_pic.py --sc 1 --rms-z 0.002 --rot-x 45.0

python sns_linac.py --sc 0
python sns_linac.py --sc 1 --dist kv
python sns_linac.py --sc 1 --dist waterbag
python sns_linac.py --sc 1 --dist gauss
python sns_linac_speed.py

python sns_ring.py --sc 0
python sns_ring.py --sc 1
python sns_ring.py --sc 1 --tilt 45.0
python sns_ring.py --sc 1 --sol 1
python sns_ring_speed.py --sc 0
python sns_ring_speed.py --sc 1
python sns_ring_orbit_correction.py
