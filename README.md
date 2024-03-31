# PyOrbit3 package installation with meson

This is experimental setup that is work in progress. There segfaulting issues with RF nodes.

## 0. Required software

One needs compilers and python.


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
Edit meson.build and set correct paths/flags for python/fftw3 headers and libraries

## 2. Build

After you have installed everything, the next step is to build. In order to build the project, navigate to the root pyorbit directory and run the following:

```bash
 pip install --no-build-isolation --editable .
```

You need only build the project after a change is made to the core c++ or python classes.

## 3. Run examples

Navigate to your examples/meson directory:

```bash
cd examples/meson
python imports_test.py
python uspas_test.py
```
