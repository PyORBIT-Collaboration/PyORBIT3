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
    "_orbit",
    sources=src,
    include_dirs=include,
    extra_compile_args=["-DUSE_MPI=0"],
)

packages = ["orbit"]
for folder in os.walk("py/orbit"):
    path = os.path.normpath(folder[0])
    path = path.split(os.sep)
    packages.append("py" + ".".join(path[1:]))

# Define the setup parameters
setup(
    ext_modules=[extension_mod],
    package_dir={
        "orbit": "src/libmain/orbit",
        "pyorbit": "py/orbit",
    },
    packages=packages,
)