import orbit


def test_import_orbit_exposes_all_submodules():
    expected_submodules = {
        "aperture",
        "bumps",
        "bunch_generators",
        "bunch_utils",
        "collimation",
        "core",
        "diagnostics",
        "errors",
        "fieldtracker",
        "foils",
        "impedances",
        "injection",
        "kickernodes",
        "lattice",
        "matching",
        "matrix_lattice",
        "orbit_correction",
        "parsers",
        "py_linac",
        "rf_cavities",
        "space_charge",
        "teapot",
        "teapot_base",
        "time_dep",
        "utils",
    }
    assert set(orbit.__all__) == expected_submodules
    assert all(hasattr(orbit, name) for name in expected_submodules)


def test_teapot_lattice_available_from_orbit_namespace():
    lattice = orbit.teapot.TEAPOT_Lattice()
    assert isinstance(lattice, orbit.teapot.TEAPOT_Lattice)
