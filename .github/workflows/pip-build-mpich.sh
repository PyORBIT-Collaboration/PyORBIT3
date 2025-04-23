python3 -m venv .po3
sed  -i 's|PATH="\$VIRTUAL_ENV\/"bin":\$PATH"|PATH="\$VIRTUAL_ENV\/"bin":\$PATH"\nsource /etc/profile.d/modules.sh\nmodule load mpi/mpich-x86_64|' .po3/bin/activate
. .po3/bin/activate
mpirun -V
pip install -U pip
pip install -r requirements.txt
pip install -U setuptools
pip install --verbose --config-settings=setup-args="-DUSE_MPI=mpich" .