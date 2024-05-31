#include "orbit_mpi.hh"

#include "wrap_LImpedance.hh"
#include "wrap_TImpedance.hh"

static PyMethodDef impedancesMethods[] = {{NULL,NULL}};


namespace wrap_impedances
{

#ifdef __cplusplus
extern "C"
{
#endif

static struct PyModuleDef cModImped = {
    PyModuleDef_HEAD_INIT,
    "impedances", "Impedances class",
    -1,
     impedancesMethods
     };

PyMODINIT_FUNC initimpedances()
{
  //create new module
  PyObject* module = PyModule_Create(&cModImped);
  wrap_impedances::initLImpedance(module);
  wrap_impedances::initTImpedance(module);
  return module;
}

PyObject* getImpedanceType(char* name)
{
  PyObject* mod = PyImport_ImportModule("orbit.core.impedances");
  PyObject* pyType = PyObject_GetAttrString(mod, name);
  Py_DECREF(mod);
  Py_DECREF(pyType);
  return pyType;
}

#ifdef __cplusplus
}
#endif // __cplusplus

} // end of namespace wrap_impedances
