# PyOrbit3 package installation

## 0. Required software
Current version doesn't support MPI or FFT based space charge calculation.
One needs compilers and python development packages, depending on  Linux flavor the package can be called **python-dev** or **python-devel**.

This guide was tested on following configurations

| CPU           | Architecture | OS           | Python  | Compiler     |
|---------------|--------------|--------------|---------|--------------|
| Intel i7-7700 | x86_64       | RHEL 8.7     | 3.9.13  | gcc-8.5      |
|               | x86_64       | Arch         | 3.10.10 | gcc-12.2.1   |
| Apple M2      | arm64        | macOS 13.3.1 | 3.9.6   | clang-14.0.3 |


## 1. Installation

If you already have configured the ESS index url (e.g. Jupyter), you can omit the index url in the command below.

```
pip install --index-url https://artifactory.esss.lu.se/artifactory/api/pypi/pypi-virtual/simple pyorbit
```

## 3. Run simple test to make sure that import works
```python
from orbit.core import orbit_mpi, bunch, teapot_base, linac, spacecharge, orbit_utils, aperture
```
The above should give no errors. It's important to `import orbit` first.

## 4. Run SNS linac example
The python part of PyOrbit should be in PYTHONPATH

```bash
cd ../pyorbit_linac_model
python pyorbit3_sns_linac_mebt_hebt2_test.py
```
