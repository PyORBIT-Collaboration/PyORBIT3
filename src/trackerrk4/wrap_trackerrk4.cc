#include "orbit_mpi.hh"

#include "wrap_trackerrk4.hh"
#include "wrap_runge_kutta_tracker.hh"
#include "wrap_py_external_effects.hh"
#include "wrap_ext_effects_container.hh"

static PyMethodDef trackerrk4Methods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  static struct PyModuleDef cModPyDem =
  {
	PyModuleDef_HEAD_INIT,
	"trackerrk4", "Tracker RK4 classes",
	-1,
	trackerrk4Methods
  };

  PyMODINIT_FUNC inittrackerrk4(){
    //create new module
    PyObject* module = PyModule_Create(&cModPyDem);
		//add the other classes init
		wrap_trackerrk4::initRungeKuttaTracker(module);
		wrap_trackerrk4_py_external_effects::initPyExternalEffects(module);
		wrap_ext_effects_container::initExtEffectsContainer(module);

		return module;
}
  PyObject* getTrackerRK4Type(const char* name){
		PyObject* mod = PyImport_ImportModule("trackerrk4");
		PyObject* pyType = PyObject_GetAttrString(mod,name);
		Py_DECREF(mod);
		Py_DECREF(pyType);
		return pyType;
	}

#ifdef __cplusplus
}
#endif
