from setuptools import Extension, setup
from pathlib import Path
import os

# main dir is special
# we need it in include but not the actual main.cc
# libmain contains python package def, we don't want it in C++ sources

src = []

for f in Path("src").rglob("*.cc"):
    excludes = ["main/main.cc"]
    include = True
    for e in excludes:
        if str(f).endswith(e):
            include = False
    if include:
        src.append(str(f))

include = []
for folder in os.walk("src"):
    excludes = ["src", "src/libmain", "src/libmain/orbit"]
    if folder[0] not in excludes:
        include.append(folder[0])
        print(folder[0])

extension_mod = Extension(
    "orbit.core._orbit",
    sources=src,
    libraries=["fftw3"],
    include_dirs=include,
    extra_compile_args=["-DUSE_MPI=1", "-fPIC", "-lmpi", "-lmpicxx", "-Wl,--enable-new-dtags"],
    extra_link_args=["-lfftw3", "-lm", "-lmpi", "-lmpicxx", "-fPIC"],
)

packages = ["orbit.core"]
for folder in os.walk("py/orbit"):
    path = os.path.normpath(folder[0])
    path = path.split(os.sep)
    packages.append(".".join(path[1:]))

package_dir = {
    "orbit": "py/orbit",
    "orbit.core": "src/libmain/orbit",
}

# This snippet generates the package structure of the orbit.core modules
# including the __init__.py file for each module
# The purpose is to be able to load individual modules from orbit.core in a
# Pythonic fashion.
core_modules = [
    "aperture",
    "orbit_mpi",
    "trackerrk4",
    "error_base",
    "bunch",
    "teapot_base",
    "linac",
    "spacecharge",
    "orbit_utils",
    "foil",
    "collimator",
    "field_sources",
    "rfcavities",
    "impedances",
    "fieldtracker",
]
for mod in core_modules:
    packages.append(f"orbit.core.{mod}")
    package_dir.update({f"orbit.core.{mod}": "src/libmain/module_template"})

# Define the setup parameters
setup(
    ext_modules=[extension_mod],
    package_dir=package_dir,
    packages=packages,
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    scripts=["bin/pyORBIT"],
)
