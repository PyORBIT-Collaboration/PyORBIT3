#include "mpi/orbit_mpi.hh"
#include "main/pyORBIT_Object.hh"

#include "orbit/Apertures/wrap_aperture.hh"

#include "orbit/Apertures/wrap_TAperture.hh"
#include "orbit/Apertures/wrap_PhaseAperture.hh"
#include "orbit/Apertures/wrap_EnergyAperture.hh"
#include "orbit/Apertures/wrap_BaseAperture.hh"
#include "orbit/Apertures/wrap_PyBaseApertureShape.hh"
#include "orbit/Apertures/wrap_PrimitiveApertureShape.hh"
#include "orbit/Apertures/wrap_CompositeApertureShape.hh"
#include "orbit/Apertures/wrap_ConvexApertureShape.hh"

#include "orbit/wrap_bunch.hh"

#include <iostream>

namespace wrap_aperture{

#ifdef __cplusplus
extern "C" {
#endif

  static PyMethodDef ApertureModuleMethods[] = { {NULL,NULL} };

	//--------------------------------------------------
	//Initialization aperture module
	//--------------------------------------------------

  static struct PyModuleDef cModPyDem =
  {
	  PyModuleDef_HEAD_INIT,
	  "aperture", "Aperture classes collection",
	  -1,
	  ApertureModuleMethods
  };
	/* The name of the function was changed to avoid collision with PyImport magic naming */
	PyMODINIT_FUNC initaperture(){
		//create new module
		PyObject* module = PyModule_Create(&cModPyDem);
		wrap_aperture::initTAperture(module);
		wrap_phase_aperture::initPhaseAperture(module);
		wrap_energy_aperture::initEnergyAperture(module);
		wrap_base_aperture::initBaseAperture(module);
		wrap_py_base_aperture_shape::initPyBaseApertureShape(module);
		wrap_primitive_aperture_shape::initPrimitiveApertureShape(module);
		wrap_py_composite_aperture_shape::initCompositeApertureShape(module);
		wrap_convex_aperture_shape::initConvexApertureShape(module);
		return module;
	}

#ifdef __cplusplus
}
#endif


}
