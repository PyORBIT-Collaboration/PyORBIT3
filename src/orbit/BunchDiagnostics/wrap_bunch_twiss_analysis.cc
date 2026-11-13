#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch_twiss_analysis.hh"
#include "wrap_bunch.hh"

#include <iostream>
#include <limits>

#include "BunchTwissAnalysis.hh"

namespace wrap_bunch_twiss_analysis{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/**
	    Constructor for python class wrapping c++ BunchTwissAnalysis instance.
      It never will be called directly.
	*/
	static PyObject* BunchTwissAnalysis_new(PyTypeObject *type, PyObject *Py_UNUSED(args), PyObject *Py_UNUSED(kwds)){
		pyORBIT_Object* self;
		self = reinterpret_cast<pyORBIT_Object*>(type->tp_alloc(type, 0));
		if (!self)
			return nullptr;
		self->cpp_obj = nullptr;
		return reinterpret_cast<PyObject*>(self);
	}

  /** This is implementation of the __init__ method */
  static int BunchTwissAnalysis_init(pyORBIT_Object *self, PyObject *Py_UNUSED(args), PyObject *Py_UNUSED(kwds)){
		self->cpp_obj =  new BunchTwissAnalysis();
	  static_cast<BunchTwissAnalysis*>(self->cpp_obj)->setPyWrapper(
          reinterpret_cast<PyObject*>(self));
    return 0;
  }

    // Match CPython's callback erasure through a generic function pointer.
    static PyCFunction fastCallMethod(
        PyObject* (*method)(PyObject*, PyObject* const*, Py_ssize_t))
    {
        return reinterpret_cast<PyCFunction>(reinterpret_cast<void (*)(void)>(method));
    }

    static BunchTwissAnalysis* getAnalysis(PyObject* self)
    {
        return static_cast<BunchTwissAnalysis*>(
            reinterpret_cast<pyORBIT_Object*>(self)->cpp_obj);
    }

    static Bunch* getBunch(PyObject* arg)
    {
        PyObject* bunchType = wrap_orbit_bunch::getBunchType("Bunch");
        const int isBunch = PyObject_IsInstance(arg, bunchType);
        if (isBunch < 0)
            return nullptr;
        if (!isBunch)
        {
            PyErr_SetString(PyExc_TypeError, "expected a Bunch");
            return nullptr;
        }
        return static_cast<Bunch*>(reinterpret_cast<pyORBIT_Object*>(arg)->cpp_obj);
    }

    // PyLong_AsLong accepts the index protocol; check the narrower C int range too.
    static bool asInt(PyObject* arg, int& value)
    {
        const long converted = PyLong_AsLong(arg);
        if (converted == -1 && PyErr_Occurred())
            return false;
        if (converted < std::numeric_limits<int>::min() ||
            converted > std::numeric_limits<int>::max())
        {
            PyErr_SetString(PyExc_OverflowError, "Python int does not fit in a C++ int");
            return false;
        }
        value = static_cast<int>(converted);
        return true;
    }

    static PyObject* getIndexedValue(PyObject* self, PyObject* arg,
                                    double (BunchTwissAnalysis::*method)(int))
    {
        int index;
        if (!asInt(arg, index))
            return nullptr;
        return PyFloat_FromDouble((getAnalysis(self)->*method)(index));
    }

    static PyObject* BunchTwissAnalysis_analyzeBunch(PyObject* self, PyObject* arg)
    {
        Bunch* bunch = getBunch(arg);
        if (!bunch)
            return nullptr;
        getAnalysis(self)->analyzeBunch(bunch);
        Py_RETURN_NONE;
    }

    static PyObject* BunchTwissAnalysis_computeBunchMoments(
        PyObject* self, PyObject* const* args, Py_ssize_t nargs)
    {
        if (nargs != 4)
        {
            PyErr_Format(PyExc_TypeError,
                         "computeBunchMoments() takes exactly 4 arguments (%zd given)", nargs);
            return nullptr;
        }
        int order, dispersion, emitnorm;
        if (!asInt(args[1], order) || !asInt(args[2], dispersion) ||
            !asInt(args[3], emitnorm))
            return nullptr;
        Bunch* bunch = getBunch(args[0]);
        if (!bunch)
            return nullptr;
        getAnalysis(self)->computeBunchMoments(bunch, order, dispersion, emitnorm);
        Py_RETURN_NONE;
    }

    static PyObject* BunchTwissAnalysis_getCorrelation(
        PyObject* self, PyObject* const* args, Py_ssize_t nargs)
    {
        if (nargs != 2)
        {
            PyErr_Format(PyExc_TypeError,
                         "getCorrelation() takes exactly 2 arguments (%zd given)", nargs);
            return nullptr;
        }
        int i, j;
        if (!asInt(args[0], i) || !asInt(args[1], j))
            return nullptr;
        return PyFloat_FromDouble(getAnalysis(self)->getCorrelation(i, j));
    }

    static PyObject* BunchTwissAnalysis_getBunchMoment(
        PyObject* self, PyObject* const* args, Py_ssize_t nargs)
    {
        if (nargs != 2)
        {
            PyErr_Format(PyExc_TypeError,
                         "getBunchMoment() takes exactly 2 arguments (%zd given)", nargs);
            return nullptr;
        }
        int i, j;
        if (!asInt(args[0], i) || !asInt(args[1], j))
            return nullptr;
        return PyFloat_FromDouble(getAnalysis(self)->getBunchMoment(i, j));
    }

    static PyObject* BunchTwissAnalysis_getAverage(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getAverage);
    }

    static PyObject* BunchTwissAnalysis_getGlobalCount(PyObject* self, PyObject* Py_UNUSED(args))
    {
        return PyLong_FromLong(getAnalysis(self)->getGlobalCount());
    }

    static PyObject* BunchTwissAnalysis_getGlobalMacrosize(PyObject* self, PyObject* Py_UNUSED(args))
    {
        return PyFloat_FromDouble(getAnalysis(self)->getGlobalMacrosize());
    }

    static PyObject* BunchTwissAnalysis_getEmittance(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getEmittance);
    }

    static PyObject* BunchTwissAnalysis_getEmittanceNormalized(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getEmittanceNormalized);
    }

    static PyObject* BunchTwissAnalysis_getAlpha(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getAlpha);
    }

    static PyObject* BunchTwissAnalysis_getBeta(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getBeta);
    }

    static PyObject* BunchTwissAnalysis_getGamma(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getGamma);
    }

    static PyObject* BunchTwissAnalysis_getDispersion(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getDispersion);
    }

    static PyObject* BunchTwissAnalysis_getDispersionDerivative(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getDispersionDerivative);
    }

    static PyObject* BunchTwissAnalysis_getTwiss(PyObject* self, PyObject* arg)
    {
        int index;
        if (!asInt(arg, index))
            return nullptr;
        BunchTwissAnalysis* analysis = getAnalysis(self);
        return Py_BuildValue("(dddd)", analysis->getAlpha(index), analysis->getBeta(index),
                             analysis->getGamma(index), analysis->getEmittance(index));
    }

    static PyObject* BunchTwissAnalysis_getEffectiveEmittance(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getEffectiveEmittance);
    }

    static PyObject* BunchTwissAnalysis_getEffectiveAlpha(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getEffectiveAlpha);
    }

    static PyObject* BunchTwissAnalysis_getEffectiveBeta(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getEffectiveBeta);
    }

    static PyObject* BunchTwissAnalysis_getEffectiveGamma(PyObject* self, PyObject* arg)
    {
        return getIndexedValue(self, arg, &BunchTwissAnalysis::getEffectiveGamma);
    }

  //-----------------------------------------------------
  //destructor for python BunchTwissAnalysis class (__del__ method).
  //-----------------------------------------------------
  static void BunchTwissAnalysis_del(pyORBIT_Object* self){
		//std::cerr<<"The BunchTwissAnalysis __del__ has been called!"<<std::endl;
		delete static_cast<BunchTwissAnalysis*>(self->cpp_obj);
		self->ob_base.ob_type->tp_free(reinterpret_cast<PyObject*>(self));
  }

	// defenition of the methods of the python BunchTwissAnalysis wrapper class
	// they will be vailable from python level
    // clang-format off
    static PyMethodDef BunchTwissAnalysisClassMethods[] = {
        {"analyzeBunch",            BunchTwissAnalysis_analyzeBunch,                        METH_O,        "Performs the Twiss analysis of the bunch."},
        {"computeBunchMoments",     fastCallMethod(BunchTwissAnalysis_computeBunchMoments), METH_FASTCALL, "Returns the XY moments of the beam up to a prescribed order"},
        {"getCorrelation",          fastCallMethod(BunchTwissAnalysis_getCorrelation),      METH_FASTCALL, "Returns the centered correlation between coordinates (ic,jc)"},
        {"getBunchMoment",          fastCallMethod(BunchTwissAnalysis_getBunchMoment),      METH_FASTCALL, "Returns the (i,j) xy moment of the beam"},
        {"getAverage",              BunchTwissAnalysis_getAverage,                          METH_O,        "Returns the average value for coordinate with index ic"},
        {"getGlobalCount",          BunchTwissAnalysis_getGlobalCount,                      METH_NOARGS,   "Returns the total number of analysed macroparticles"},
        {"getGlobalMacrosize",      BunchTwissAnalysis_getGlobalMacrosize,                  METH_NOARGS,   "Returns the total macrosize"},
        {"getEmittance",            BunchTwissAnalysis_getEmittance,                        METH_O,        "Returns the emittance for index 0,1,2 - x,y,z planes"},
        {"getEmittanceNormalized",  BunchTwissAnalysis_getEmittanceNormalized,              METH_O,        "Returns the normalized emittance for index 0,1 - x,y planes"},
        {"getAlpha",                BunchTwissAnalysis_getAlpha,                            METH_O,        "Returns Twiss alpha for index 0,1,2 - x,y,z planes"},
        {"getBeta",                 BunchTwissAnalysis_getBeta,                             METH_O,        "Returns Twiss beta for index 0,1,2 - x,y,z planes"},
        {"getGamma",                BunchTwissAnalysis_getGamma,                            METH_O,        "Returns Twiss gamma for index 0,1,2 - x,y,z planes"},
        {"getTwiss",                BunchTwissAnalysis_getTwiss,                            METH_O,        "Returns Twiss tuple (alpha,beta,gamma,emitt) for index 0,1,2 - x,y,z planes"},
        {"getDispersion",           BunchTwissAnalysis_getDispersion,                       METH_O,        "Returns Twiss dispersion function for index 0,1 - x,y planes"},
        {"getDispersionDerivative", BunchTwissAnalysis_getDispersionDerivative,             METH_O,        "Returns Twiss dispersion' function for index 0,1 - x,y planes"},
        {"getEffectiveEmittance",   BunchTwissAnalysis_getEffectiveEmittance,               METH_O,        "Returns the effective emittance for index 0,1 - x,y planes"},
        {"getEffectiveAlpha",       BunchTwissAnalysis_getEffectiveAlpha,                   METH_O,        "Returns effective Twiss alpha for index 0,1 - x,y planes"},
        {"getEffectiveBeta",        BunchTwissAnalysis_getEffectiveBeta,                    METH_O,        "Returns effective Twiss beta for index 0,1 - x,y planes"},
        {"getEffectiveGamma",       BunchTwissAnalysis_getEffectiveGamma,                   METH_O,        "Returns effective Twiss gamma for index 0,1 - x,y planes"},
        {nullptr, nullptr, 0, nullptr}
    };
    // clang-format on

	// defenition of the memebers of the python BunchTwissAnalysis wrapper class
	// they will be vailable from python level
	static PyMemberDef BunchTwissAnalysisClassMembers [] = {
		{NULL}
	};

	//new python BunchTwissAnalysis wrapper type definition
	static PyTypeObject pyORBIT_BunchTwissAnalysis_Type = {
		PyVarObject_HEAD_INIT(NULL, 0)
		"BunchTwissAnalysis", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		reinterpret_cast<destructor>(BunchTwissAnalysis_del) , /*tp_dealloc*/
		0, /*tp_print*/
		0, /*tp_getattr*/
		0, /*tp_setattr*/
		0, /*tp_compare*/
		0, /*tp_repr*/
		0, /*tp_as_number*/
		0, /*tp_as_sequence*/
		0, /*tp_as_mapping*/
		0, /*tp_hash */
		0, /*tp_call*/
		0, /*tp_str*/
		0, /*tp_getattro*/
		0, /*tp_setattro*/
		0, /*tp_as_buffer*/
		Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
		"The BunchTwissAnalysis python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		BunchTwissAnalysisClassMethods, /* tp_methods */
		BunchTwissAnalysisClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		reinterpret_cast<initproc>(BunchTwissAnalysis_init), /* tp_init */
		0, /* tp_alloc */
		BunchTwissAnalysis_new, /* tp_new */
	};



	//--------------------------------------------------
	//Initialization BunchTwissAnalysis of the pyBunchTwissAnalysis class
	//--------------------------------------------------
  void initbunchtwissanalysis(PyObject* module){
		if (PyType_Ready(&pyORBIT_BunchTwissAnalysis_Type) < 0) return;
		Py_INCREF(&pyORBIT_BunchTwissAnalysis_Type);
		PyModule_AddObject(module, "BunchTwissAnalysis", reinterpret_cast<PyObject*>(&pyORBIT_BunchTwissAnalysis_Type));
	}

#ifdef __cplusplus
}
#endif


}
