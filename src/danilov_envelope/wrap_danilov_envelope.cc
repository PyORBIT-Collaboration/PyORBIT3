#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_DanilovEnvelopeSolver20.hh"
#include "wrap_DanilovEnvelopeSolver22.hh"
#include "wrap_bunch.hh"
#include "wrap_danilov_envelope.hh"

namespace wrap_danilov_envelope {

#ifdef __cplusplus
extern "C" {
#endif

static PyMethodDef DanilovEnvelopeModuleMethods[] = {{NULL, NULL}};

static struct PyModuleDef cModPyDem = {
  PyModuleDef_HEAD_INIT, 
  "danilov_envelope", "Danilov distribution envelope solvers", 
  -1, 
  DanilovEnvelopeModuleMethods
};

// The name of the function was changed to avoid collision with PyImport magic naming.
PyMODINIT_FUNC initdanilovenvelope() {
  // create new module
  PyObject *module = PyModule_Create(&cModPyDem);
  wrap_danilov_envelope::initDanilovEnvelopeSolver20(module);
  wrap_danilov_envelope::initDanilovEnvelopeSolver22(module);
  return module;
}

#ifdef __cplusplus
}
#endif

} // namespace wrap_danilov_envelope