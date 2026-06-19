import pytest
from orbit.core.bunch import Bunch


@pytest.fixture
def bunch():
    """Create a basic proton bunch."""
    b = Bunch()
    b.mass(0.938272)
    b.charge(1.0)
    b.getSyncParticle().kinEnergy(1.0)
    return b


class TestBunchCreation:

    def test_empty_bunch_has_no_particles(self, bunch):
        assert bunch.getSize() == 0
        assert bunch.getTotalCount() == 0

    def test_default_mass_and_charge(self, bunch):
        assert bunch.mass() == pytest.approx(0.938272)
        assert bunch.charge() == pytest.approx(1.0)


class TestAddParticle:

    def test_add_particle_stores_coords_correctly(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.x(idx) == 1.0
        assert bunch.px(idx) == 0.01
        assert bunch.y(idx) == 2.0
        assert bunch.py(idx) == 0.02
        assert bunch.z(idx) == 3.0
        assert bunch.pz(idx) == 0.03

    def test_add_multiple_particles_increases_size(self, bunch):
        for i in range(5):
            bunch.addParticle(
                float(i) * 1e-3, float(i) * 1e-4,
                float(i) * 1e-3, float(i) * 1e-4,
                float(i) * 1e-3, float(i) * 1e-4,
            )
        assert bunch.getSize() == 5

    def test_add_particle_returns_consecutive_indices(self, bunch):
        indices = []
        for i in range(5):
            idx = bunch.addParticle(
                float(i) * 1e-3, float(i) * 1e-4,
                float(i) * 1e-3, float(i) * 1e-4,
                float(i) * 1e-3, float(i) * 1e-4,
            )
            indices.append(idx)
        assert indices == [0, 1, 2, 3, 4]

    def test_many_particles_dont_interfere(self, bunch):
        """Each particle's coordinates are independent."""
        for i in range(10):
            bunch.addParticle(float(i), 0.0, 0.0, 0.0, 0.0, 0.0)
        for i in range(10):
            assert bunch.x(i) == float(i)


class TestCoordinateAccessors:

    def test_get_x(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.x(idx) == 1.0

    def test_get_px(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.px(idx) == 0.01

    def test_get_y(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.y(idx) == 2.0

    def test_get_py(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.py(idx) == 0.02

    def test_get_z(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.z(idx) == 3.0

    def test_get_pz(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.pz(idx) == 0.03

    def test_set_x(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.x(idx, 99.0)
        assert bunch.x(idx) == 99.0

    def test_set_y(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.y(idx, 88.0)
        assert bunch.y(idx) == 88.0

    def test_set_z(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.z(idx, 77.0)
        assert bunch.z(idx) == 77.0

    def test_set_px(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.px(idx, 0.99)
        assert bunch.px(idx) == 0.99

    def test_set_py(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.py(idx, 0.88)
        assert bunch.py(idx) == 0.88

    def test_set_pz(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.pz(idx, 0.77)
        assert bunch.pz(idx) == 0.77

    def test_xp_is_px(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.xp(idx) == bunch.px(idx)

    def test_yp_is_py(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.yp(idx) == bunch.py(idx)

    def test_dE_is_pz(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        assert bunch.dE(idx) == bunch.pz(idx)

    def test_set_xp(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.xp(idx, 0.99)
        assert bunch.xp(idx) == 0.99
        assert bunch.px(idx) == 0.99

    def test_set_yp(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.yp(idx, 0.88)
        assert bunch.yp(idx) == 0.88
        assert bunch.py(idx) == 0.88

    def test_set_dE(self, bunch):
        idx = bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.dE(idx, 0.77)
        assert bunch.dE(idx) == 0.77
        assert bunch.pz(idx) == 0.77


class TestParticleLifecycle:

    def test_delete_particle(self, bunch):
        for i in range(5):
            bunch.addParticle(float(i) * 1e-3, 0, 0, 0, 0, 0)
        assert bunch.getSize() == 5
        bunch.deleteParticle(2)
        assert bunch.getSize() == 4

    def test_delete_particle_fast_no_compress(self, bunch):
        for i in range(5):
            bunch.addParticle(float(i) * 1e-3, 0, 0, 0, 0, 0)
        bunch.deleteParticleFast(2)
        assert bunch.getSize() == 5
        assert bunch.flag(2) == 0

    def test_delete_particle_fast_then_compress(self, bunch):
        for i in range(5):
            bunch.addParticle(float(i) * 1e-3, 0, 0, 0, 0, 0)
        bunch.deleteParticleFast(2)
        bunch.compress()
        assert bunch.getSize() == 4

    def test_flag_for_alive(self, bunch):
        idx = bunch.addParticle(0, 0, 0, 0, 0, 0)
        assert bunch.flag(idx) == 1

    def test_flag_for_dead(self, bunch):
        idx = bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.deleteParticleFast(idx)
        assert bunch.flag(idx) == 0

    def test_recover_particle(self, bunch):
        idx = bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.deleteParticleFast(idx)
        assert bunch.flag(idx) == 0
        bunch.recoverParticle(idx)
        assert bunch.flag(idx) == 1

    def test_delete_all_particles(self, bunch):
        for i in range(10):
            bunch.addParticle(0, 0, 0, 0, 0, 0)
        assert bunch.getSize() == 10
        bunch.deleteAllParticles()
        assert bunch.getSize() == 0
        assert bunch.getTotalCount() == 0

    def test_compress_no_dead_particles(self, bunch):
        for i in range(3):
            bunch.addParticle(float(i), 0, 0, 0, 0, 0)
        size_before = bunch.getSize()
        bunch.compress()
        assert bunch.getSize() == size_before
        for i in range(3):
            assert bunch.x(i) == float(i)

    def test_compress_multiple_times(self, bunch):
        for i in range(5):
            bunch.addParticle(float(i) * 1e-3, 0, 0, 0, 0, 0)
        bunch.deleteParticleFast(2)
        bunch.compress()
        size1 = bunch.getSize()
        bunch.compress()
        assert bunch.getSize() == size1

    def test_hard_delete_then_recover(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.addParticle(1.0, 0, 0, 0, 0, 0)
        assert bunch.getSize() == 2
        bunch.deleteParticle(0)
        assert bunch.getSize() == 1
        assert bunch.x(0) == 1.0
        bunch.recoverParticle(0)
        assert bunch.flag(0) == 1

    def test_flag_after_hard_delete(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.deleteParticle(0)
        assert bunch.flag(0) == 1

    def test_delete_last_particle(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.deleteParticle(0)
        assert bunch.getSize() == 0


class TestBunchAttributes:

    def test_mass(self, bunch):
        assert bunch.mass() == pytest.approx(0.938272)
        bunch.mass(0.139570)
        assert bunch.mass() == pytest.approx(0.139570)

    def test_charge(self, bunch):
        assert bunch.charge() == pytest.approx(1.0)
        bunch.charge(-1.0)
        assert bunch.charge() == pytest.approx(-1.0)

    def test_macro_size(self, bunch):
        bunch.macroSize(1e10)
        assert bunch.macroSize() == 1e10

    def test_classical_radius(self, bunch):
        rp = bunch.classicalRadius()
        assert rp == pytest.approx(1.534698e-18, rel=1e-6)

    def test_b_rho_positive(self, bunch):
        assert bunch.B_Rho() > 0

    def test_bunch_attr_double_round_trip(self, bunch):
        bunch.bunchAttrDouble("my_double", 42.0)
        assert bunch.bunchAttrDouble("my_double") == 42.0

    def test_bunch_attr_int_round_trip(self, bunch):
        bunch.bunchAttrInt("my_int", 42)
        assert bunch.bunchAttrInt("my_int") == 42

    def test_has_bunch_attr_double(self, bunch):
        assert bunch.hasBunchAttrDouble("mass") == 1
        assert bunch.hasBunchAttrDouble("nonexistent") == 0

    def test_bunch_attr_names(self, bunch):
        names = bunch.bunchAttrDoubleNames()
        assert "mass" in names
        assert "charge" in names
        assert "classical_radius" in names
        assert "macro_size" in names

    def test_macro_size_default(self, bunch):
        assert bunch.macroSize() == 0.0

    def test_bunch_attr_int_names(self, bunch):
        bunch.bunchAttrInt("my_int", 42)
        names = bunch.bunchAttrIntNames()
        assert "my_int" in names

    def test_has_bunch_attr_int(self, bunch):
        assert bunch.hasBunchAttrInt("nonexistent") == 0
        bunch.bunchAttrInt("my_int", 42)
        assert bunch.hasBunchAttrInt("my_int") == 1

    def test_b_rho_after_energy_change(self, bunch):
        brho1 = bunch.B_Rho()
        sp = bunch.getSyncParticle()
        sp.kinEnergy(2.0)
        assert bunch.B_Rho() != pytest.approx(brho1, rel=1e-6)

    def test_classical_radius_after_mass_change(self, bunch):
        rc1 = bunch.classicalRadius()
        bunch.mass(0.139570)
        rc2 = bunch.classicalRadius()
        assert rc2 != rc1


class TestParticleAttributes:

    def test_add_part_attr(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.addPartAttr("ParticleIdNumber")
        assert bunch.hasPartAttr("ParticleIdNumber") == 1

    def test_part_attr_value_round_trip(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.addPartAttr("ParticleIdNumber")
        bunch.partAttrValue("ParticleIdNumber", 0, 0, 42)
        assert bunch.partAttrValue("ParticleIdNumber", 0, 0) == 42

    def test_remove_part_attr(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.addPartAttr("ParticleIdNumber")
        assert bunch.hasPartAttr("ParticleIdNumber") == 1
        bunch.removePartAttr("ParticleIdNumber")
        assert bunch.hasPartAttr("ParticleIdNumber") == 0

    def test_get_part_attr_names(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.addPartAttr("ParticleIdNumber")
        names = bunch.getPartAttrNames()
        assert "ParticleIdNumber" in names

    def test_remove_all_part_attr(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.addPartAttr("TurnNumber")
        bunch.addPartAttr("ParticleIdNumber")
        assert len(bunch.getPartAttrNames()) == 2
        bunch.removeAllPartAttr()
        assert bunch.getPartAttrNames() == ()

    def test_get_part_attr_dicts(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.addPartAttr("TurnNumber")
        dicts = bunch.getPartAttrDicts()
        assert "TurnNumber" in dicts
        assert "size" in dicts["TurnNumber"]

    def test_get_possible_part_attr_names(self, bunch):
        names = bunch.getPossiblePartAttrNames()
        assert "TurnNumber" in names
        assert "ParticleIdNumber" in names

    def test_get_part_attr_size(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.addPartAttr("TurnNumber")
        assert bunch.getPartAttrSize("TurnNumber") == 1

    def test_multiple_particles_with_attributes(self, bunch):
        bunch.addPartAttr("ParticleIdNumber")
        for i in range(5):
            bunch.addParticle(float(i) * 1e-3, 0, 0, 0, 0, 0)
            bunch.partAttrValue("ParticleIdNumber", i, 0, 100 + i)
        for i in range(5):
            assert bunch.partAttrValue(
                "ParticleIdNumber", i, 0
            ) == pytest.approx(100.0 + i)


class TestCopyOperations:

    def test_copy_empty_bunch_to(self, bunch):
        src = bunch
        src.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        dst = Bunch()
        src.copyEmptyBunchTo(dst)
        assert dst.getSize() == 0
        assert dst.mass() == pytest.approx(src.mass())
        assert dst.charge() == pytest.approx(src.charge())

    def test_copy_bunch_to(self, bunch):
        src = bunch
        src.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        dst = Bunch()
        src.copyBunchTo(dst)
        assert dst.getSize() == 1
        assert dst.x(0) == 1.0
        assert dst.y(0) == 2.0
        assert dst.z(0) == 3.0

    def test_add_particles_to(self, bunch):
        src = bunch
        src.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        dst = Bunch()
        dst.mass(0.938272)
        dst.charge(1.0)
        src.addParticlesTo(dst)
        assert dst.getSize() == 1
        assert dst.x(0) == 1.0
        assert dst.py(0) == 0.02


class TestRingWrap:

    def test_ringwrap_inside_stays(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 1.0, 0)
        bunch.ringwrap(20.0)
        assert bunch.z(0) == pytest.approx(1.0)

    def test_ringwrap_wraps(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 15.0, 0)
        bunch.ringwrap(20.0)
        assert bunch.z(0) == pytest.approx(-5.0)

    def test_ringwrap_negative(self, bunch):
        bunch.addParticle(0, 0, 0, 0, -15.0, 0)
        bunch.ringwrap(20.0)
        assert bunch.z(0) == pytest.approx(5.0)

    def test_ringwrap_at_zero(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 0, 0)
        bunch.ringwrap(20.0)
        assert bunch.z(0) == pytest.approx(0.0)

    def test_ringwrap_beyond_full_ring(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 25.0, 0)
        bunch.ringwrap(20.0)
        assert bunch.z(0) == pytest.approx(-5.0)
        bunch.addParticle(0, 0, 0, 0, -25.0, 0)
        bunch.ringwrap(20.0)
        assert bunch.z(1) == pytest.approx(5.0)

    def test_ringwrap_idempotent(self, bunch):
        bunch.addParticle(0, 0, 0, 0, 15.0, 0)
        bunch.ringwrap(20.0)
        z_first = bunch.z(0)
        bunch.ringwrap(20.0)
        assert bunch.z(0) == pytest.approx(z_first)


class TestBunchCapacity:

    def test_get_capacity(self, bunch):
        assert bunch.getCapacity() > 0

    def test_get_size_global(self, bunch):
        assert bunch.getSizeGlobal() == 0

    def test_get_size_global_from_memory(self, bunch):
        assert bunch.getSizeGlobalFromMemory() == 0

    def test_get_total_count_tracks_allocation(self, bunch):
        for i in range(5):
            bunch.addParticle(0, 0, 0, 0, 0, 0)
        assert bunch.getTotalCount() == 5
        bunch.deleteParticleFast(2)
        bunch.addParticle(99, 0, 0, 0, 0, 0)
        assert bunch.getTotalCount() == 6


class TestSyncParticle:

    def test_get_sync_particle(self, bunch):
        sp = bunch.getSyncParticle()
        assert sp.mass() == pytest.approx(0.938272)

    def test_sync_particle_kin_energy(self, bunch):
        sp = bunch.getSyncParticle()
        sp.kinEnergy(2.0)
        assert sp.kinEnergy() == pytest.approx(2.0)

    def test_sync_particle_beta_gamma(self, bunch):
        sp = bunch.getSyncParticle()
        beta = sp.beta()
        gamma = sp.gamma()
        assert gamma == pytest.approx(1.0 / (1.0 - beta**2) ** 0.5, rel=1e-6)
        assert beta > 0
        assert gamma > 1.0


class TestFileIO:

    def test_dump_and_read_bunch_empty(self, bunch, tmp_path):
        fname = str(tmp_path / "empty.bunch")
        bunch.dumpBunch(fname)
        b2 = Bunch()
        b2.mass(0.938272)
        b2.charge(1.0)
        b2.getSyncParticle().kinEnergy(1.0)
        b2.readBunch(fname)
        assert b2.getSize() == 0
        assert b2.mass() == pytest.approx(0.938272)
        assert b2.charge() == pytest.approx(1.0)

    def test_dump_and_read_bunch_single(self, bunch, tmp_path):
        fname = str(tmp_path / "single.bunch")
        bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.dumpBunch(fname)
        b2 = Bunch()
        b2.mass(0.938272)
        b2.charge(1.0)
        b2.getSyncParticle().kinEnergy(1.0)
        b2.readBunch(fname)
        assert b2.getSize() == 1
        assert b2.x(0) == 1.0
        assert b2.px(0) == 0.01
        assert b2.y(0) == 2.0
        assert b2.py(0) == 0.02
        assert b2.z(0) == 3.0
        assert b2.pz(0) == 0.03

    def test_dump_and_read_bunch_multi(self, bunch, tmp_path):
        fname = str(tmp_path / "multi.bunch")
        for i in range(3):
            bunch.addParticle(float(i), float(i) * 0.01, float(i) * 2, float(i) * 0.02, float(i) * 3, float(i) * 0.03)
        bunch.dumpBunch(fname)
        b2 = Bunch()
        b2.mass(0.938272)
        b2.charge(1.0)
        b2.getSyncParticle().kinEnergy(1.0)
        b2.readBunch(fname)
        assert b2.getSize() == 3
        for i in range(3):
            assert b2.x(i) == float(i)
            assert b2.px(i) == float(i) * 0.01
            assert b2.y(i) == float(i) * 2
            assert b2.pz(i) == float(i) * 0.03

    def test_dump_and_read_with_particle_attrs(self, bunch, tmp_path):
        fname = str(tmp_path / "attrs.bunch")
        bunch.addParticle(1.0, 0.01, 2.0, 0.02, 3.0, 0.03)
        bunch.addPartAttr("TurnNumber")
        bunch.partAttrValue("TurnNumber", 0, 0, 42)
        bunch.dumpBunch(fname)
        b2 = Bunch()
        b2.mass(0.938272)
        b2.charge(1.0)
        b2.getSyncParticle().kinEnergy(1.0)
        b2.readBunch(fname)
        assert b2.getSize() == 1
        assert b2.x(0) == 1.0
        assert b2.hasPartAttr("TurnNumber") == 1
        assert b2.partAttrValue("TurnNumber", 0, 0) == 42
