from setuptools import Extension, setup
from pathlib import Path


def list_cc(path):
    return [str(f) for f in Path(path).iterdir() if f.is_file() and f.name.endswith('.cc')]


def list_dirs(path, exclude=None):
    exclude_list = [] if exclude is None else exclude
    p = Path(path)
    return [p] + [d for d in p.rglob('*')if d.is_dir() and d not in exclude_list]


top_level_dirs = ['mpi', 'utils', 'orbit', 'linac', 'spacecharge', 'teapot']
all_dirs = [d for sublist in [list_dirs(f'src/{p}') for p in top_level_dirs] for d in sublist]
all_src = [list_cc(f'{p}') for p in all_dirs]

# main dir is special
# we need it in include but not the actual main.cc
# libmain contains python package def, we don't want it in C++ sources
include = [f'{d}' for d in all_dirs] + ['src/main']
src = [d for sublist in all_src for d in sublist] +['src/libmain/libmain.cc']


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
                   # 'orbit': 'py/orbit',
                   },
      packages=['orbit3',
                # 'orbit',
                ],
      )

