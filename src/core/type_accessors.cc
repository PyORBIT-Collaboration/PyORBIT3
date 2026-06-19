#include <Python.h>

#include "wrap_bunch.hh"
#include "wrap_spacecharge.hh"
#include "wrap_mpi_comm.hh"
#include "wrap_utils.hh"
#include "wrap_trackerrk4.hh"

extern "C" {

namespace wrap_orbit_bunch {

PyObject* getBunchType(const char* name){
    PyObject* mod = PyImport_ImportModule("orbit.core.bunch");
    PyObject* pyType = PyObject_GetAttrString(mod,name);
    Py_DECREF(mod);
	Py_DECREF(pyType);
    return pyType;
}

}

PyObject* getSpaceChargeType(const char* name){
    PyObject* mod = PyImport_ImportModule("orbit.core.spacecharge");
    PyObject* pyType = PyObject_GetAttrString(mod,name);
    Py_DECREF(mod);
    Py_DECREF(pyType);
    return pyType;
}

namespace wrap_orbit_mpi_comm {

PyObject* getMPI_CommType(const char* name){
    PyObject* mod = PyImport_ImportModule("orbit.core.orbit_mpi");
    PyObject* mpi_comm_mod = PyObject_GetAttrString(mod,"mpi_comm");
    PyObject* pyType = PyObject_GetAttrString(mpi_comm_mod,name);
    Py_DECREF(mpi_comm_mod);
    Py_DECREF(mod);
    Py_DECREF(pyType);
    return pyType;
}

}

namespace wrap_orbit_utils {

PyObject* getOrbitUtilsType(const char* name){
    PyObject* mod = PyImport_ImportModule(const_cast<char*>("orbit.core.orbit_utils"));
    PyObject* pyType = PyObject_GetAttrString(mod,name);
    Py_DECREF(mod);
    Py_DECREF(pyType);
    return pyType;
}

}

PyObject* getTrackerRK4Type(const char* name){
    PyObject* mod = PyImport_ImportModule("orbit.core.trackerrk4");
    PyObject* pyType = PyObject_GetAttrString(mod,name);
    Py_DECREF(mod);
    Py_DECREF(pyType);
    return pyType;
}

}
