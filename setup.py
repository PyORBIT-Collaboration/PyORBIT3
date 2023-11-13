from setuptools import Extension, setup
import os
from pathlib import Path
import subprocess


def parse_options(line):
    library_dirs = []
    libraries = []
    include_dirs = []
    compile_options = []
    for option in line.split(' ')[1:]:
        if option.startswith('-L'):
            library_dirs.append(option[2:])
        elif option.startswith('-l'):
            libraries.append(option[2:])
        elif option.startswith('-I'):
            include_dirs.append(option[2:])
        else:
            compile_options.append(option[:])
    return library_dirs, libraries, include_dirs, compile_options


library_dirs, libraries, include_dirs, extra_compile_args = [], [], [], []


mpi_compiler = os.environ.get('MPICC', default='mpicc')

if mpi_compiler:
    try:
        result = subprocess.check_output([mpi_compiler, '-showme']).decode().strip()
        options = parse_options(result)
        library_dirs, libraries, include_dirs, extra_compile_args = options
    except Exception:
        print("MPICC not valid")
        library_dirs, libraries, include_dirs, extra_compile_args = [], [], [], ["-DUSE_MPI=0"]
        print('MPI not found will be compiled without MPI support')

fftw_include = os.environ.get('FFTW3_INCLUDE_DIR', default=None)
if fftw_include:
    include_dirs.append(fftw_include)
fftw_lib = os.environ.get('FFTW3_LIB_DIR', default=None)
if fftw_lib:
    library_dirs.append(fftw_lib)
libraries += ['fftw3']


src = []

for f in Path("src").rglob("*.cc"):
    excludes = ["main/main.cc"]
    include = True
    for e in excludes:
        if str(f).endswith(e):
            include = False
    if include:
        src.append(str(f))


for folder in os.walk("src"):
    excludes = ["src", "src/libmain", "src/libmain/orbit"]
    if folder[0] not in excludes:
        include_dirs.append(folder[0])
        print(folder[0])

extension_mod = Extension(
    "orbit.core._orbit",
    sources=src,
    libraries=libraries,
    include_dirs=include_dirs,
    library_dirs=library_dirs,
    # unclear when next line is needed
    runtime_library_dirs=library_dirs,
    extra_compile_args=extra_compile_args + ["-Wl,--enable-new-dtags"],
    extra_link_args=["-lm", "-fPIC"],
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

