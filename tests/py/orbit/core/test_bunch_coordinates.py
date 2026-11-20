import pytest

from orbit.core.bunch import Bunch

np = pytest.importorskip("numpy")
pytestmark = pytest.mark.skipif(
    not hasattr(Bunch, "coordinates"), reason="requires USE_NUMPY=true"
)


@pytest.mark.parametrize("count", [0, 1, 1000])
def test_coordinates_shape_values_and_round_trip(count):
    expected = np.arange(count * 6, dtype=np.float64).reshape(count, 6)
    bunch = Bunch()
    for row in expected:
        bunch.addParticle(*row)

    actual = bunch.coordinates()
    assert actual.shape == (count, 6)
    assert actual.dtype == np.float64
    assert actual.flags.c_contiguous
    assert actual.flags.owndata
    np.testing.assert_array_equal(actual, expected)

    other = Bunch()
    other.addParticles(actual)
    np.testing.assert_array_equal(other.coordinates(), expected)


def test_coordinates_is_independent_snapshot():
    bunch = Bunch()
    bunch.addParticle(*range(6))
    snapshot = bunch.coordinates()
    snapshot[0, 0] = -100
    assert bunch.x(0) == 0
    bunch.x(0, 42)
    assert snapshot[0, 0] == -100
    bunch.deleteAllParticles()
    np.testing.assert_array_equal(snapshot, [[-100, 1, 2, 3, 4, 5]])
