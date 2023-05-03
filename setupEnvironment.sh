command_exists () {
    type "$1" &> /dev/null ;
}


export ORBIT3_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "ORBIT installed in $ORBIT3_ROOT"

export ORBIT3_ARCH=`uname -s`

# alias python='python2'

if command_exists python3; then
   PYEX=python3
else
   PYEX=python
fi

export PYTHON3_VERSION=`$PYEX -c "from distutils import sysconfig; print(sysconfig.get_config_var('VERSION'));"`
echo "Python version is $PYTHON3_VERSION"

PYTHON3_LIB_DIR=`$PYEX -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBPL'));"`
#if [ -f $PYTHON3_LIB_DIR/libpython${PYTHON3_VERSION}.a ]
#then
#  export PYTHON3_ROOT_LIB=$PYTHONi3_LIB_DIR/libpython${PYTHON3_VERSION}.a
#  LIB_TYPE=static
#else
  export PYTHON3_ROOT_LIB="-L $PYTHON3_LIB_DIR -lpython${PYTHON3_VERSION}"
  LIB_TYPE=dynamic
#fi

echo "Found python library: ${PYTHON3_LIB_DIR} will use $LIB_TYPE library"

export PYTHON3_ROOT_INC=`$PYEX -c "from distutils import sysconfig; print(sysconfig.get_config_var('INCLUDEPY'));"`
echo "Found Python include directory: $PYTHON3_ROOT_INC"

export PYTHONPATH=${PYTHONPATH}:${ORBIT_ROOT}/py:${ORBIT_ROOT}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ORBIT_ROOT}/lib

if command_exists mpirun ; then
   echo "Found mpirun at: `which mpirun`"
   MPI_RUN_DIR=`dirname $(which mpirun)`
else
    MPI_RUN_DIR=`dirname $(find /usr 2>/dev/null| fgrep bin/mpirun | head -n1)`
    export PATH=$PATH:$MPI_RUN_DIR
    echo "Added  $MPI_RUN_DIR to PATH"
fi

export MPI_CPP=$MPI_RUN_DIR/mpicxx
echo "MPI_CPP set to $MPI_CPP"
