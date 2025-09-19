#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_danilov_20_envelope_solver.hh"
#include "wrap_danilov_22_envelope_solver.hh"
#include "wrap_bunch.hh"
#include "wrap_envelope.hh"

namespace wrap_envelope {

#ifdef __cplusplus
extern "C" {
#endif

static PyMethodDef EnvelopeModuleMethods[] = {{NULL, NULL}};

static struct PyModuleDef cModPyDem = {
  PyModuleDef_HEAD_INIT, 
  "envelope", "Beam envelope solvers", 
  -1, 
  EnvelopeModuleMethods
};

// Name changed to avoid collision with PyImport magic naming.
PyMODINIT_FUNC initenvelope() {
  PyObject *module = PyModule_Create(&cModPyDem);
  wrap_envelope::initdanilov20envelopesolver(module);
  wrap_envelope::initdanilov22envelopesolver(module);
  return module;
}

#ifdef __cplusplus
}
#endif

} // namespace wrap_envelope