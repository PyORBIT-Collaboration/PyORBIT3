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
    excludes = ["src", "src/libmain", "src/libmain/orbit3"]
    if folder[0] not in excludes:
        include.append(folder[0])
        print(folder[0])

extension_mod = Extension('_orbit3',
                          sources=src,
                          include_dirs=include,
                          extra_compile_args=['-DUSE_MPI=0'],
                          )

# Define the setup parameters
setup(name='orbit3',
      version='1.0',
      description='A C++ extension module for Python.',
      ext_modules=[extension_mod],
      package_dir={'orbit3': 'src/libmain/orbit3',
                   'orbit': 'py/orbit',
                   },
      packages=['orbit3',
                'orbit',
                ],
      )

