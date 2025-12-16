Installation
=============================

1. Introduction
---------------

This guide provides instructions how to install PyORBIT code. This guide
doesn’t cover MPI enabled installation. The following configurations are
included in CI testing, versions will change as the runner images
progress.

+------------+-----------+------------+-------+-----------+-----------+
| HW         | Arc       | OS         | P     | Compiler  | Package   |
|            | hitecture |            | ython |           |           |
+============+===========+============+=======+===========+===========+
| PC         | x86_64    | CentOS     | 3     | g         | pip-24.0  |
|            |           | latest     | .9.18 | cc-11.4.1 |           |
+------------+-----------+------------+-------+-----------+-----------+
| PC         | x86_64    | Ubuntu     | 3     | g         | pip-24.0  |
|            |           | latest     | .12.3 | cc-13.2.0 |           |
+------------+-----------+------------+-------+-----------+-----------+
| Apple      | arm64     | macOS 14   | 3     | cla       | pip-24.0  |
| Silicon    |           |            | .12.3 | ng-15.0.0 |           |
+------------+-----------+------------+-------+-----------+-----------+
| PC         | x86_64    | Ubuntu     | 3.    | g         | con       |
|            |           | latest     | 10.14 | cc-13.2.0 | da-24.5.0 |
+------------+-----------+------------+-------+-----------+-----------+

2. Installation from source
---------------------------

First step is to clone the source code:

.. code:: bash

   git clone https://github.com/PyORBIT-Collaboration/PyORBIT3.git

Pip Setup
~~~~~~~~~

**PIP** based setup is more involved, we recommend using **conda** if
unsure. #### Prepare OS The following commands may require root access.

.. raw:: html

   <details>

.. raw:: html

   <summary>

Click for Ubuntu-based distributions

.. raw:: html

   </summary>

.. code:: bash

     apt-get update
     apt-get install -y  build-essential python3 libfftw3-dev python3-venv libpython3-dev pkg-config git

.. raw:: html

   </details>

.. raw:: html

   <details>

.. raw:: html

   <summary>

Click for Redhat-based distributions

.. raw:: html

   </summary>

.. code:: bash

     dnf group install -y "Development Tools"
     dnf install -y python3-devel fftw3-devel

.. raw:: html

   </details>

.. raw:: html

   <details>

.. raw:: html

   <summary>

Click for Mac

.. raw:: html

   </summary>

Install Homebrew, make sure that homebrew programs are in the **$PATH**
(optional step in Homebrew installation)
``bash   brew install pkg-config fftw``

.. raw:: html

   </details>

#### Create Python virtual environment Make sure that you have the
correct python version installed. We require python>3.9. Create virtual
environment.
``python3 -m venv .po3   . .po3/bin/activate   pip install -U pip   pip install -r requirements.txt   pip install -U setuptools``

Conda Setup
~~~~~~~~~~~

First of all make sure you have conda installed and development
packages. Development packages for Ubuntu:

::

   apt update -y
   apt install -y curl gpg git build-essential

Then run the following:

.. code:: bash

   cd pyorbit3
   conda env create -n po3 --file environment.yml
   conda activate po3
   pip install -U meson-python setuptools setuptools-scm

3. Build
--------

If you plan to modify PyORBIT’s code, install it in editable mode. You
will NOT need to rebuild after modifications to the code.
`Meson <MesonBuild.md>`__ will rebuild as necessary on import.

::

   pip install --no-build-isolation --editable .

Alternatively if you don’t plan to modify PyORBIT’s code

::

   pip install .

4. Run full SNS linac example
-----------------------------

Navigate to your **examples** directory and launch tracking of SNS
linac.

.. code:: bash

   cd examples/SNS_Linac/pyorbit3_linac_model/
   python pyorbit3_sns_linac_mebt_hebt2.py

5. MPI consideration
--------------------

By default, the build system will try to find MPI and compile against
it. You can control which MPI to use with command line option when
building.

.. code:: bash

   pip install --config-settings=setup-args="-DUSE_MPI=none" .

Above will build PyORBIT without MPI even if MPI is present. You can
change that option to ``mpich``, ``ompi``, ``none`` or ``auto``
(default). \| MPI flavor \| Installation command \| \|—————\|————–\| \|
No MPI \|
``pip install --config-settings=setup-args="-DUSE_MPI=none" .`` \| \|
The first found MPI installation if any \|
``pip install --config-settings=setup-args="-DUSE_MPI=auto" .`` \| \|
OpenMPI \|
``pip install --config-settings=setup-args="-DUSE_MPI=ompi" .`` \| \|
MPICH \|
``pip install --config-settings=setup-args="-DUSE_MPI=mpich" .`` \|

Meson uses PKG_CONFIG to discover packages. It could be useful to help
it to find your MPI installation:

.. code:: bash

   PKG_CONFIG_PATH=/opt/lib/pkgconfig pip install --verbose .
