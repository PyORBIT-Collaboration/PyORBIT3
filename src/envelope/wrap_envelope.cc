#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch.hh"
#include "wrap_envelope.hh"
#include "wrap_env_solver_kv.hh"

namespace wrap_envelope {

#ifdef __cplusplus
extern "C" {
#endif

static PyMethodDef envelopeModuleMethods[] = {{NULL, NULL}};

static struct PyModuleDef cModPyDem = {
  PyModuleDef_HEAD_INIT, 
  "envelope", "Beam envelope solvers", 
  -1, 
  envelopeModuleMethods
};

PyMODINIT_FUNC initenvelope() {
  PyObject *module = PyModule_Create(&cModPyDem);
  wrap_envelope::initEnvSolverKV(module);
  return module;
}

#ifdef __cplusplus
}
#endif

}
