from pathlib import Path

import pytest

from orbit.teapot import TEAPOT_Lattice


def test_sns_ring_teapot():
    madx_file = Path(__file__).parent / "inputs" / "sns_ring.lat"

    lattice = TEAPOT_Lattice()
    lattice.readMADX(str(madx_file), "rnginjsol")
    lattice.initialize()

    lattice_length_madx = 248.0098418
    assert abs(lattice_length_madx - lattice.getLength()) < 1e-7
