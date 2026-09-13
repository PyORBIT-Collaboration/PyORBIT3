#include "mpi/orbit_mpi.hh"

#include "linac/wrap_linacmodule.hh"
#include "linac/rfgap/wrap_BaseRfGap.hh"
#include "linac/rfgap/wrap_BaseRfGap_slow.hh"
#include "linac/rfgap/wrap_MatrixRfGap.hh"
#include "linac/rfgap/wrap_RfGapTTF.hh"
#include "linac/rfgap/wrap_RfGapTTF_slow.hh"
#include "linac/rfgap/wrap_SuperFishFieldSource.hh"
#include "linac/rfgap/wrap_RfGapThreePointTTF.hh"
#include "linac/rfgap/wrap_RfGapThreePointTTF_slow.hh"
#include "linac/tracking/wrap_linac_tracking.hh"

static PyMethodDef linacmoduleMethods[] = { {NULL,NULL} };

static struct PyModuleDef linacModDef =
        {
                PyModuleDef_HEAD_INIT,
                "linac", "Linac C++ classes",
                -1,
                linacmoduleMethods
        };


#ifdef __cplusplus
extern "C" {
#endif
  namespace wrap_linac{

    PyMODINIT_FUNC initlinac(){
        //create new module
        PyObject* module = PyModule_Create(&linacModDef);
        //add the other classes init
        wrap_linac::initBaseRfGap(module);
        wrap_linac::initBaseRfGap_slow(module);
        wrap_linac::initMatrixRfGap(module);
        wrap_linac::initRfGapTTF(module);
        wrap_linac::initRfGapTTF_slow(module);
        wrap_linac::initSuperFishFieldSource(module);
        wrap_linac::initRfGapThreePointTTF(module);
        wrap_linac::initRfGapThreePointTTF_slow(module);
        //initialization of the linac tracking module
        wrap_linac_tracking::initlinactracking(module);
        return module;
    }

    PyObject* getLinacType(char* name){
        PyObject* mod = PyImport_ImportModule("orbit.core.linac");
        PyObject* pyType = PyObject_GetAttrString(mod,name);
        Py_DECREF(mod);
        Py_DECREF(pyType);
        return pyType;
    }
  }

#ifdef __cplusplus
}
#endif
