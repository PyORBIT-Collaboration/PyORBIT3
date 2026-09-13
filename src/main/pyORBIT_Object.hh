//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    pyORBIT_Object.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    06/08/2007
//
// DESCRIPTION
//    The wrapper around any pointer to the ORBIT C++ class instance.
//
///////////////////////////////////////////////////////////////////////////

#ifndef PY_ORBIT_OBJECT_H
#define PY_ORBIT_OBJECT_H

#include <Python.h>
#include <structmember.h>

#ifdef __cplusplus
extern "C" {
#endif

 typedef struct {
   PyObject_HEAD
   void* cpp_obj;
 } pyORBIT_Object;

#ifdef __cplusplus
}

namespace pyorbit {

// Registry entries are borrowed references owned by their Python wrappers.
void registerPyWrapper(void* cpp_obj, PyObject* py_wrapper);
void unregisterPyWrapper(void* cpp_obj, PyObject* py_wrapper);
PyObject* getPyWrapper(void* cpp_obj);

}
#endif

#endif
