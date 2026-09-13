#include "mpi/orbit_mpi.hh"

#include "utils/wrap_utils.hh"
#include "utils/matrix/wrap_matrix.hh"
#include "utils/matrix/wrap_phase_vector.hh"
#include "utils/wrap_py_base_field_source.hh"
#include "utils/wrap_field_source_container.hh"
#include "utils/wrap_function.hh"
#include "utils/wrap_splinech.hh"
#include "utils/statistics/wrap_statmoments2d.hh"
#include "utils/bunch/wrap_bunch_extrema_calculator.hh"
#include "utils/integration/wrap_gauss_legendre_integrator.hh"
#include "utils/polynomial/wrap_polynomial.hh"
#include "utils/wrap_numrecipes.hh"
#include "utils/bunch/wrap_bunch_utils_functions.hh"
#include "utils/harmonic_analysis/wrap_harmonic_data.hh"
#include "utils/wrap_random.hh"
//#include "wrap_field_sources_module.hh"

namespace wrap_orbit_utils{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

  static PyMethodDef UtilsModuleMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

	static struct PyModuleDef cModPyDem =
	{
		PyModuleDef_HEAD_INIT,
		"orbit_utils", "PyORBIT utilities.",
		-1,
		UtilsModuleMethods
	};

  PyMODINIT_FUNC initutils(){
    //create new module
    PyObject* module = PyModule_Create(&cModPyDem);
		//add the other classes init
		wrap_utils_martix::initMatrix(module);
		wrap_utils_phase_vector::initPhaseVector(module);
		wrap_utils_py_base_field_source::initPyBaseFieldSource(module);
		wrap_field_source_container::initFieldSourceContainer(module);
		wrap_function::initFunction(module);
		wrap_splinech::initSplineCH(module);
		wrap_statmoments2d::initstatmoments2d(module);
		wrap_utils_bunch::initBunchExtremaCalculator(module);
		wrap_gl_integrator::initGLIntegrator(module);
		wrap_polynomial::initPolynomial(module);
		//this call creates the module orbit_utils.numrecipes with functions
		wrap_numrecipes::initNumrecipes(module);
		wrap_utils_bunch_functions::initBunchUtilsFunctions(module);
		wrap_harmonicdata::initHarmonicData(module);
		//wrap_field_sources_module::initFieldSourcesModule(module);
		wrap_random::initRandom(module);
		return module;
  }

#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_utils
}
