#ifndef PY_BASE_APERTURE_SHAPE_H
#define PY_BASE_APERTURE_SHAPE_H

#include <Python.h>

#include "orbit/Bunch.hh"
#include "orbit/Apertures/BaseApertureShape.hh"


///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   PyBaseApertureShape
//
// AUTHOR:
//   Andrei Shishlo October 2022
//
///////////////////////////////////////////////////////////////////////////

/**
  PyBaseApertureShape is an implementation of BaseApertureShape class which
  will be exposed to Python level. It uses the Python methods to calculate
  the result of int PyBaseApertureShape::inside(...) method, so it is very
  slow and should be used only for testing, developing prototypes etc.tc.
*/

class PyBaseApertureShape: public BaseApertureShape
{
public:

	/** PyBaseApertureShape constructor */
	PyBaseApertureShape(PyObject* py_wrapper);

	/** PyBaseApertureShape decstructor */
	virtual ~PyBaseApertureShape();

	/** Return 1 if the particular macro-particle is inside this shape */
	int inside(Bunch* bunch, int count);

private:
	PyObject* py_wrapper;

};

//end of PY_BASE_APERTURE_SHAPE_H ifdef
#endif
