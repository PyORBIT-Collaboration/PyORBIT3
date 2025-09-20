#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch.hh"
#include "wrap_envelope.hh"
#include "wrap_kv_envelope_tracker.hh"
#include "wrap_danilov_envelope_tracker.hh"


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
  wrap_envelope::initKVEnvelopeTracker(module);
  wrap_envelope::initDanilovEnvelopeTracker(module);
  return module;
}

#ifdef __cplusplus
}
#endif

}
