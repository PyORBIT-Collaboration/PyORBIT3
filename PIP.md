## Alternative PIP Setup
You don't need it if you used conda.

PIP setup requires several libraries installed with your package manager.
The package names are different for different OS flavors/types.
Here are required packages for Debian based (tested on the latest Ubuntu), Redhat based (tested on CentOS Stream 9) and Mac (tested on Apple Silicon).


On Debian based distributions:
```
apt-get update
apt-get install build-essential python3 openmpi-bin openmpi-common libopenmpi-dev
apt-get install zlib1g-dev libfftw3-dev python3-distutils python3-venv libpython3-dev git
```

On RedHat based distributions
```
dnf group install "Development Tools"
dnf install python3-devel openmpi-devel fftw3-devel
```

On Mac
You need to have [HomeBrew](https://brew.sh/) installed.

```
brew install fftw
```

Prepare virtual environment with *install/create_env.y* script.

```
python3 install/create_env.py <PATH_TO_YOUR_ENV>
```
It makes sense to put your virtual environment outside PyORBIT directory.
The script should find MPI and FFTW and add them to corresponding environment variables. These variables are set only when your virtual enviroment is activated.

Currently, it will need an MPI distribution that supports **mpicc -showme** option. 
MPICH doesn't support it, OPENMPI does. So if you have MPICH installed along with OPENMPI you need to specify that you will need to use OPENMPI.

```
python3 install/create_env.py --mpi /usr/lib64/openmpi <PATH_TO_YOUR_ENV>
```

Activate your environment
```
. <PATH_TO_YOUR_ENV>/bin/activate
```

You are ready to [build](README.md#2-build)
