# PyOrbit3 package installation

## 1. Introduction
This guide provides instructions how to install PyORBIT code. <br>
This guide doesn't cover MPI enabled installation. <br>
The following configurations are included in CI testing, versions will change as the runner images progress.

| HW            | Architecture | OS            | Python  | Compiler     | Package      |
|---------------|--------------|---------------|---------|--------------|--------------|
| PC            | x86_64       | CentOS latest | 3.9.18  | gcc-11.4.1   | pip-24.0     |
| PC            | x86_64       | Ubuntu latest | 3.12.3  | gcc-13.2.0   | pip-24.0     |
| Apple Silicon | arm64        | macOS 14      | 3.12.3  | clang-15.0.0 | pip-24.0     |
| PC            | x86_64       | Ubuntu latest | 3.10.14 | gcc-13.2.0   | conda-24.5.0 |



## 2. Installation from source

First step is to clone the source code:

```bash
git clone https://github.com/PyORBIT-Collaboration/PyORBIT3.git
```

### Pip Setup

#### Ubuntu based distributions:
```
sudo apt-get update
sudo apt-get install -y  build-essential python3 libfftw3-dev python3-venv libpython3-dev pkg-config git
```

#### Redhat based distributions:
```
dnf group install -y "Development Tools"
dnf install -y python3-devel fftw3-devel
```

#### MacOS
Install Homebrew, make sure that  homebrew programs are in the **$PATH** (optional step in Homebrew installation)
```bash
brew install pkg-config fftw
```

Make sure that you have the correct python version installed. We require python>3.9. <br>
Create virtual environment.
```
python3 -m venv .po3
. .po3/bin/activate
pip install -U pip
pip install -r requirements.txt
pip install -U setuptools
```

### Conda Setup

First of all make sure you have conda installed and development packages.<br>
Development packages for Ubuntu:
```
apt update -y
apt install -y curl gpg git build-essential
```

Then run the following:

```bash
cd pyorbit3
conda env create -n po3 --file environment.yml
conda activate po3
pip install -U meson-python setuptools setuptools-scm
```


## 3. Build

If you plan to modify PyORBIT's code, install it in editable mode. 
You will NOT need to rebuild after modifications to the code. [Meson](MesonBuild.md) will rebuild as necessary on import.
```
pip install --no-build-isolation --editable .
```

Alternatively if you don't plan to modify PyORBIT's code
```
pip install .
```


## 4. Run full SNS linac example

Navigate to your **examples** directory and launch tracking of SNS linac.

```bash
cd examples/SNS_Linac/pyorbit3_linac_model/
python pyorbit3_sns_linac_mebt_hebt2.py
```



