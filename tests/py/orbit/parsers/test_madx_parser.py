from pathlib import Path

import pytest

from orbit.parsers import MADX_Parser
from orbit.teapot import TEAPOT_Lattice


def test_repeated_sequence_element_creates_independent_occurrences(tmp_path):
    file = Path(__file__).parent / "inputs" / "repeated_element.lat"

    parser = MADX_Parser()
    parser.parse(str(file))
    sequence = parser.getSequenceList()

    bpm_elements = [element for element in sequence if element.getName() == "bpm"]
    assert len(bpm_elements) == 2
    assert bpm_elements[0] is not bpm_elements[1]
    assert [element.getParameter("position") for element in bpm_elements] == [1.0, 2.0]

    probe = next(element for element in sequence if element.getName() == "probe")
    assert probe.getParameter("position") == pytest.approx(2.5)

    drift_lengths = [
        element.getParameter("l")
        for element in sequence
        if element.getType() == "drift"
    ]
    assert drift_lengths == pytest.approx([1.0, 1.0, 0.5, 0.5])
    assert sum(element.getParameter("l") for element in sequence) == pytest.approx(3.0)


def test_sns_ring_teapot():
    madx_file = Path(__file__).parent / "inputs" / "sns_ring.lat"

    lattice = TEAPOT_Lattice()
    lattice.readMADX(str(madx_file), "rnginjsol")
    lattice.initialize()

    lattice_length_madx = 248.0098418
    assert abs(lattice_length_madx - lattice.getLength()) < 1e-7
