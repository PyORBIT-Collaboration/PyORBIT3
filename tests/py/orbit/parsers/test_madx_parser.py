import pytest

from orbit.parsers import MADX_Parser


def test_repeated_sequence_element_creates_independent_occurrences(tmp_path):
    madx_file = tmp_path / "repeated_element.madx"
    madx_file.write_text(
        """bpm: monitor;
probe: marker;
ring: sequence, l=3.0;
bpm, at=1.0;
bpm, at=2.0;
probe, at=0.5, from=bpm;
endsequence;
""",
        encoding="utf-8",
    )

    parser = MADX_Parser()
    parser.parse(str(madx_file))
    sequence = parser.getSequenceList()

    bpm_occurrences = [element for element in sequence if element.getName() == "bpm"]
    assert len(bpm_occurrences) == 2
    assert bpm_occurrences[0] is not bpm_occurrences[1]
    assert [element.getParameter("position") for element in bpm_occurrences] == [1.0, 2.0]

    probe = next(element for element in sequence if element.getName() == "probe")
    assert probe.getParameter("position") == pytest.approx(2.5)

    drift_lengths = [
        element.getParameter("l")
        for element in sequence
        if element.getType() == "drift"
    ]
    assert drift_lengths == pytest.approx([1.0, 1.0, 0.5, 0.5])
    assert sum(element.getParameter("l") for element in sequence) == pytest.approx(3.0)
