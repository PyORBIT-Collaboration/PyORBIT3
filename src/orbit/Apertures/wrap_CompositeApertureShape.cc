#include "mpi/orbit_mpi.hh"
#include "main/pyORBIT_Object.hh"

#include "orbit/wrap_bunch.hh"

#include <iostream>

#include "orbit/Apertures/CompositeApertureShape.hh"

namespace wrap_py_composite_aperture_shape{
	typedef struct {
		PyObject_HEAD
		void* cpp_obj;
		PyObject* shapes;
	} pyORBIT_CompositeApertureShape;

#ifdef __cplusplus
extern "C" {
#endif

	/**
	Constructor for python class wrapping c++ Circle, Ellipse, and Rectangular ApertureShape instances.
	*/
	static PyObject* CompositeApertureShape_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_CompositeApertureShape* self;
		self = (pyORBIT_CompositeApertureShape *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		self->shapes = NULL;
		return (PyObject *) self;
	}

  /** This is implementation of the __init__ method */
  static int CompositeApertureShape_init(pyORBIT_CompositeApertureShape *self, PyObject *args, PyObject *kwds){
  	self->cpp_obj = (pyORBIT_Object*) new CompositeApertureShape();
	  self->shapes = PyList_New(0);
	  if(self->shapes == NULL) {
		  delete ((CompositeApertureShape*) self->cpp_obj);
		  self->cpp_obj = NULL;
		  return -1;
	  }
	  pyorbit::registerPyWrapper(self->cpp_obj, (PyObject*) self);
    return 0;
  }

	// name([name]) - sets or returns the name of the shape
  static PyObject* CompositeApertureShape_name(PyObject *self, PyObject *args){
    pyORBIT_Object* pyCompositeApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyCompositeApertureShape->cpp_obj;
    const char* name = NULL;
    if(!PyArg_ParseTuple(	args,"|s:name",&name)){
      ORBIT_MPI_Finalize("CompositeApertureShape.name(...) - call should be - name([name]). Stop.");
    }
		if(name != NULL){
      std::string name_str(name);
      cpp_BaseApertureShape->setName(name_str);
		}
		return Py_BuildValue("s",cpp_BaseApertureShape->getName().c_str());
  }

	// typeName() - returns the type name of the shape
  static PyObject* CompositeApertureShape_typeName(PyObject *self, PyObject *args){
    pyORBIT_Object* pyCompositeApertureShape= (pyORBIT_Object*) self;
		BaseApertureShape* cpp_BaseApertureShape = (BaseApertureShape*) pyCompositeApertureShape->cpp_obj;
		int nVars = PyTuple_Size(args);
		if(nVars != 0){
			ORBIT_MPI_Finalize("CompositeApertureShape.typeName() - has no parameters. Stop.");
    }
		return Py_BuildValue("s",cpp_BaseApertureShape->getTypeName().c_str());
  }

	// addApertureShape() - adds the ApertureShape instance to composite
  static PyObject* CompositeApertureShape_addApertureShape(PyObject *self, PyObject *args){
		pyORBIT_CompositeApertureShape* composite = (pyORBIT_CompositeApertureShape*) self;
		CompositeApertureShape* cpp_CompositeApertureShape = (CompositeApertureShape*) composite->cpp_obj;
	  PyObject* pyBaseApertureShape;
		if(!PyArg_ParseTuple(args,"O:setApertureShape",&pyBaseApertureShape)){
				ORBIT_MPI_Finalize("CompositeApertureShape.addApertureShape(BaseApertureShape) - parameter is needed. Stop.");
		}
		if(PyList_Append(composite->shapes, pyBaseApertureShape) < 0) return NULL;
		cpp_CompositeApertureShape->addApertureShape((BaseApertureShape*) ((pyORBIT_Object*) pyBaseApertureShape)->cpp_obj);
		Py_INCREF(Py_None);
		return Py_None;
  }

	// getApertureShapes() - returns the ApertureShape instances inside the composite
  static PyObject* CompositeApertureShape_getApertureShapes(PyObject *self, PyObject *args){
		pyORBIT_CompositeApertureShape* composite = (pyORBIT_CompositeApertureShape*) self;
		return PyList_AsTuple(composite->shapes);
  }

  //-----------------------------------------------------
  //destructor for python CompositeApertureShape class (__del__ method).
  //-----------------------------------------------------
  static void CompositeApertureShape_del(pyORBIT_CompositeApertureShape* self){
		//std::cerr<<"debug CompositeApertureShape __del__ has been called!"<<std::endl;
		pyorbit::unregisterPyWrapper(self->cpp_obj, (PyObject*) self);
		delete ((BaseApertureShape*)self->cpp_obj);
		Py_CLEAR(self->shapes);
		self->ob_base.ob_type->tp_free((PyObject*)self);
  }

	// definition of the methods of the python CompositeApertureShape wrapper class
	// they will be vailable from python level
	static PyMethodDef CompositeApertureShapeClassMethods[] = {
		{ "name",                 CompositeApertureShape_name,            METH_VARARGS,"Sets or returns the name of the shape."},
		{ "typeName",             CompositeApertureShape_typeName,        METH_VARARGS,"Returns the type of the shape."},
		{ "addApertureShape",     CompositeApertureShape_addApertureShape,METH_VARARGS,"Adds the ApertureShape instance to composit."},
		{ "getApertureShapes", CompositeApertureShape_getApertureShapes,  METH_VARARGS,"Returns tuple of ApertureShape instances."},
		{NULL}
  };

	static PyMemberDef CompositeApertureShapeClassMembers [] = {
		{NULL}
	};


	//new python CompositeApertureShape wrapper type definition
	static PyTypeObject pyORBIT_CompositeApertureShape_Type = {
		PyVarObject_HEAD_INIT(NULL, 0)
		"CompositeApertureShape", /*tp_name*/
		sizeof(pyORBIT_CompositeApertureShape), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) CompositeApertureShape_del , /*tp_dealloc*/
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
		"The CompositeApertureShape python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		CompositeApertureShapeClassMethods, /* tp_methods */
		CompositeApertureShapeClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) CompositeApertureShape_init, /* tp_init */
		0, /* tp_alloc */
		CompositeApertureShape_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization CompositeApertureShape class
	//--------------------------------------------------

	void initCompositeApertureShape(PyObject* module){
		//check that the CompositeApertureShape wrapper is ready
		if (PyType_Ready(&pyORBIT_CompositeApertureShape_Type) < 0) return;
		Py_INCREF(&pyORBIT_CompositeApertureShape_Type);
		PyModule_AddObject(module, "CompositeApertureShape", (PyObject *)&pyORBIT_CompositeApertureShape_Type);
	}

#ifdef __cplusplus
}
#endif


}
