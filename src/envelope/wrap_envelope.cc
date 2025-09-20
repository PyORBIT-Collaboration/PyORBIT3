#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch.hh"
#include "wrap_envelope.hh"
#include "wrap_env_solver_danilov_20.hh"
#include "wrap_env_solver_danilov_22.hh"

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

// The name of the function was changed to avoid collision with PyImport magic naming.
PyMODINIT_FUNC initenvelope() {
  // create new module
  PyObject *module = PyModule_Create(&cModPyDem);
  wrap_envelope::initEnvSolverDanilov20(module);
  wrap_envelope::initEnvSolverDanilov22(module);
  return module;
}

#ifdef __cplusplus
}
#endif

} // namespace wrap_envelope