#include <Python.h>

#include "orbit/wrap_bunch.hh"
#include "spacecharge/wrap_spacecharge.hh"
#include "mpi/wrap_mpi_comm.hh"
#include "utils/wrap_utils.hh"
#include "trackerrk4/wrap_trackerrk4.hh"

extern "C" {

namespace {

PyObject* getModuleAttribute(const char* module_name, const char* name){
    PyObject* module = PyImport_ImportModule(module_name);
    if(module == NULL){
        return NULL;
    }
    PyObject* attribute = PyDict_GetItemString(PyModule_GetDict(module), name);
    Py_DECREF(module);
    if(attribute == NULL){
        PyErr_Format(PyExc_AttributeError, "module '%s' has no attribute '%s'", module_name, name);
    }
    return attribute;
}

}

namespace wrap_orbit_bunch {

PyObject* getBunchType(const char* name){
    return getModuleAttribute("orbit.core.bunch", name);
}

}

PyObject* getSpaceChargeType(const char* name){
    return getModuleAttribute("orbit.core.spacecharge", name);
}

namespace wrap_orbit_mpi_comm {

PyObject* getMPI_CommType(const char* name){
    PyObject* mpi_comm_module = getModuleAttribute("orbit.core.orbit_mpi", "mpi_comm");
    if(mpi_comm_module == NULL){
        return NULL;
    }
    PyObject* type = PyDict_GetItemString(PyModule_GetDict(mpi_comm_module), name);
    if(type == NULL){
        PyErr_Format(PyExc_AttributeError, "module 'mpi_comm' has no attribute '%s'", name);
    }
    return type;
}

}

namespace wrap_orbit_utils {

PyObject* getOrbitUtilsType(const char* name){
    return getModuleAttribute("orbit.core.orbit_utils", name);
}

}

PyObject* getTrackerRK4Type(const char* name){
    return getModuleAttribute("orbit.core.trackerrk4", name);
}

}
