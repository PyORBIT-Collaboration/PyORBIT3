#include "orbit_mpi.hh"

#include "wrap_danilov_20_envelope_solver.hh"
#include "wrap_danilov_22_envelope_solver.hh"

static PyMethodDef envelopeMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

static struct PyModuleDef cModPyDem =
{
	PyModuleDef_HEAD_INIT,
	"envelope", "Beam envelope solver classes",
	-1,
	envelopeMethods
};

PyMODINIT_FUNC initenvelope(){
	PyObject* module = PyModule_Create(&cModPyDem);
	wrap_envelope::initDanilov20EnvelopeSolver(module);
	wrap_envelope::initDanilov22EnvelopeSolver(module);
	return module;
}

PyObject* getenvelopeType(const char* name){
	PyObject* mod = PyImport_ImportModule("orbit.core.envelope");
	PyObject* pyType = PyObject_GetAttrString(mod,name);
	Py_DECREF(mod);
	Py_DECREF(pyType);
	return pyType;
}

#ifdef __cplusplus
}
#endif