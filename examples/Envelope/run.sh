#!/bin/bash

set -x 

python test_env_2d_fodo.py
python test_env_2d_fodo_speed.py
python test_env_3d_drift.py

cd sns_linac
python test_sns_linac.py
cd ..

cd sns_ring
python test_sns_ring.py
python test_sns_ring_speed.py
cd ..
