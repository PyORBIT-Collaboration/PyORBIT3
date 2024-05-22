# PyOrbit3 with meson

This uses meson-python to build orbit package.

There is no **setup.py** file, instead we have **meson.build**.
**pyproject.toml** is changed to use meson.

This is experimental setup that is work in progress.
The pure python part is built with hierarchical **meson.build** files in **py/**.
The C++ setup is combined in one file **src/meson.build**.

### Main modifications in C++ code
1. **src/libmain/** is not used, still there for reference but will be gone soon.
2. **src/core/** contains one C++ file per module inside _orbit.core_
3. The files **wrap_XXXX.cc** were modified to correctly reference modules 
```cpp
// line
PyObject* mod = PyImport_ImportModule("_bunch");
// replaced with 
PyObject* mod = PyImport_ImportModule("orbit.core.bunch");
```



# Setup

## 0. Required software

One needs compilers, python and libfftw (and potentially mpi).
See [PyORBIT3](https://github.com/PyORBIT-Collaboration/PyORBIT3) for external 
requirements. 


## 1. Preparing environment

First step is to clone the source code from meson branch:

```bash
git clone -b meson https://github.com/azukov/PyORBIT3.git
```

Initialize new virtual environment and install packages

```
python -m venv .mes
source .mes/bin/activate
pip install -U pip 
pip install -r requirements
```
Edit **meson.build** and set correct paths/flags for python/fftw3 headers and libraries

## 2. Build

To install orbit package in development mode run following:
```bash
 pip install --no-build-isolation --editable .
```
No rebuild is necessary, just edit **py/** or **src/** and meson will rebuild as needed when import happens.


## 3. Run examples

Special examples used for meson testing

```bash
cd examples/meson
python imports_test.py
python uspas_test.py
```

SNS linac example
```bash
cd examples/SNS_Linac/pyorbit3_linac_model/
python pyorbit3_sns_linac_mebt_hebt2.py
```
