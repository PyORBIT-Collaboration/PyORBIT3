#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"

#include <iostream>
#include <string>

#include "Random.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_random
{

	void error(const char *msg) { ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
	extern "C"
	{
#endif

		static PyObject *random_seed(PyObject *self, PyObject *args)
		{
			int n;
			if (!PyArg_ParseTuple(args, "i:seed_n", &n))
			{
				error("seed(n) - parameters are needed.");
			}
			Random::seed(n);

			return Py_None;
		}

		static PyObject *random_ran1(PyObject *self, PyObject *args)
		{
			return Py_BuildValue("d", Random::ran1());
		}

		static PyMethodDef RandomModuleMethods[] = {
			{"seed", random_seed, METH_VARARGS, "seed(n) seed the random number generator."},
			{"ran1", random_ran1, METH_VARARGS, "ran1() - generate a random number uniformly distributed between 0 and 1."},
			{NULL, NULL, 0, NULL} /* Sentinel */
		};

		static struct PyModuleDef cModPyDem =
			{
				PyModuleDef_HEAD_INIT,
				"random", "Random number generator using the C++11 std library.",
				-1,
				RandomModuleMethods};

		//--------------------------------------------------
		// Initialization functions of the random module
		//--------------------------------------------------
		void initRandom(PyObject *module)
		{
			// create random module
			PyObject *module_ran = PyModule_Create(&cModPyDem);
			Py_INCREF(module_ran);
			PyModule_AddObject(module, const_cast<char *>("random"), module_ran);
		}

#ifdef __cplusplus
	}
#endif

}
