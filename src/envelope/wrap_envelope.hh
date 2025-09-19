#include "orbit_mpi.hh"

#include "wrap_danilov_20_envelope_solver.hh"
#include "wrap_danilov_21_envelope_solver.hh"
#include "wrap_envelope.hh"
#include "wrap_envelopecalc2p5d.hh"
#include "wrap_envelopeforcecalc2p5d.hh"
#include "wrap_envelopecalc2p5d_rb.hh"
#include "wrap_envelopecalc_slicebyslice_2D.hh"
#include "wrap_lspacechargecalc.hh"
#include "wrap_envelopecalc3d.hh"
#include "wrap_uniform_ellipsoid_field_calculator.hh"
#include "wrap_envelopecalc_uniform_ellipse.hh"

static PyMethodDef spacechargeMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  static struct PyModuleDef cModPyDem =
  {
	  PyModuleDef_HEAD_INIT,
	  "spacecharge", "Space Charge classes",
	  -1,
	  spacechargeMethods
  };

  PyMODINIT_FUNC initspacecharge(){
    //create new module
    PyObject* module = PyModule_Create(&cModPyDem);
		//add the other classes init
		wrap_envelope::initGrid1D(module);
		wrap_envelope::initGrid2D(module);
		wrap_envelope::initGrid3D(module);
		wrap_envelope::initUniformEllipsoidFieldCalculator(module);
		wrap_envelope::initSpaceChargeCalcUniformEllipse(module);
		wrap_envelope::initPoissonSolverFFT2D(module);
		wrap_envelope::initPoissonSolverFFT3D(module);
		wrap_envelope::initBoundary2D(module);
		wrap_envelope::initSpaceChargeCalc2p5D(module);
		wrap_envelope::initSpaceChargeCalc2p5Drb(module);
		wrap_envelope::initSpaceChargeCalcSliceBySlice2D(module);
		wrap_lspacechargecalc::initLSpaceChargeCalc(module);
		wrap_envelope::initSpaceChargeCalc3D(module);
		wrap_envelope::initSpaceChargeForceCalc2p5D(module);
		return module;
  }

	PyObject* getEnvelopeType(const char* name){
		PyObject* mod = PyImport_ImportModule("orbit.core.envelope");
		PyObject* pyType = PyObject_GetAttrString(mod,name);
		Py_DECREF(mod);
		Py_DECREF(pyType);
		return pyType;
	}

#ifdef __cplusplus
}
#endif
