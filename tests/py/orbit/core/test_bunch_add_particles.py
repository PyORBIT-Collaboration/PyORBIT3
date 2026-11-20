import pytest

from orbit.core.bunch import Bunch

np = pytest.importorskip("numpy")
pytestmark = pytest.mark.skipif(
    not hasattr(Bunch, "addParticles"), reason="requires USE_NUMPY=true"
)


def coordinates(bunch):
    getters = (bunch.x, bunch.px, bunch.y, bunch.py, bunch.z, bunch.dE)
    return np.array([[get(i) for get in getters] for i in range(bunch.getSize())])


@pytest.mark.parametrize("data", [
    np.arange(18).reshape(3, 6).tolist(),
    tuple(tuple(range(i, i + 6)) for i in range(3)),
    np.arange(18, dtype=np.float64).reshape(3, 6),
    np.arange(18, dtype=np.float32).reshape(3, 6),
    np.arange(18, dtype=np.int64).reshape(3, 6),
    np.arange(36, dtype=float).reshape(3, 12)[:, ::2],
    np.asfortranarray(np.arange(18.).reshape(3, 6)),
    np.arange(18.).reshape(3, 6)[::-1],
    np.arange(18.).reshape(3, 6).astype('>f8'),
])
def test_add_particles_matches_scalar_append(data):
    actual, expected = Bunch(), Bunch()
    for bunch in (actual, expected):
        bunch.addParticle(*range(6))
    assert actual.addParticles(data) is None
    for row in data:
        expected.addParticle(*row)
    np.testing.assert_array_equal(coordinates(actual), coordinates(expected))


def test_empty_array_is_noop():
    bunch = Bunch()
    bunch.addParticle(*range(6))
    assert bunch.addParticles(np.empty((0, 6))) is None
    np.testing.assert_array_equal(coordinates(bunch), [range(6)])


@pytest.mark.parametrize("data", [
    [], list(range(6)), np.zeros((6, 2)), np.zeros((2, 6, 1)), 1,
    [[0] * 6, [0] * 5], [[0] * 6, ["bad"] * 6],
])
def test_invalid_input_does_not_modify_bunch(data):
    bunch = Bunch()
    bunch.addParticle(*range(6))
    with pytest.raises((TypeError, ValueError)):
        bunch.addParticles(data)
    np.testing.assert_array_equal(coordinates(bunch), [range(6)])


def test_wrong_argument_count():
    bunch = Bunch()
    with pytest.raises(TypeError):
        bunch.addParticles()
    with pytest.raises(TypeError):
        bunch.addParticles([], [])


def test_particle_attributes_and_compression_match_scalar_append():
    actual, expected = Bunch(), Bunch()
    for bunch in (actual, expected):
        bunch.addPartAttr("ParticleIdNumber")
        bunch.addParticle(*range(6))
        bunch.partAttrValue("ParticleIdNumber", 0, 0, 42)
        bunch.addParticle(*range(6))
        bunch.deleteParticleFast(1)
    data = np.arange(12.).reshape(2, 6)
    actual.addParticles(data)
    for row in data:
        expected.addParticle(*row)
    assert actual.getSize() == expected.getSize()
    for bunch in (actual, expected):
        bunch.compress()
    np.testing.assert_array_equal(coordinates(actual), coordinates(expected))
    for i in range(expected.getSize()):
        assert actual.partAttrValue("ParticleIdNumber", i, 0) == expected.partAttrValue(
            "ParticleIdNumber", i, 0
        )
    assert actual.partAttrValue("ParticleIdNumber", 0, 0) == 42
