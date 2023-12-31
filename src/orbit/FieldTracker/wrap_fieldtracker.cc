#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#include "wrap_fieldtracker.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "FieldTracker.hh"

#include "wrap_fieldtracker.hh"

namespace wrap_fieldtracker{

    void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
    extern "C" {
#endif

        /**
         Constructor for python class wrapping c++ FieldTracker instance.
         It never will be called directly.
         */
        static PyObject* FieldTracker_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
            pyORBIT_Object* self;
            self = (pyORBIT_Object *) type->tp_alloc(type, 0);
            self->cpp_obj = NULL;
            return (PyObject *) self;
        }

        /** This is implementation of the __init__ method */
        static int FieldTracker_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
            double bx = 0.0;
            double by = 0.0;
            double ax = 0.0;
            double ay = 0.0;
            double ex = 0.0;
            double epx = 0.0;
            double l = 0.0;
            double zi = 0.0;
            double zf = 0.0;
            double ds = 0.0;
            int niters = 0;
            double resid = 0.0;
            double xrefi = 0.0;
            double yrefi = 0.0;
            double eulerai = 0.0;
            double eulerbi = 0.0;
            double eulergi = 0.0;
            FieldTracker* cpp_FieldTracker = (FieldTracker*)((pyORBIT_Object*) self)->cpp_obj;
            PyObject* pyBunch;

            const char* filename = NULL;

            //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
            if(!PyArg_ParseTuple( args,"ddddddddddiddddddOs:arguments", &bx, &by,
         	       &ax, &ay, &ex,  &epx, &l, &zi,  &zf, &ds,  &niters, &resid,
         	       &xrefi, &yrefi, &eulerai, &eulerbi, &eulergi, &pyBunch, &filename))
            {
                error("PyBunch - addParticle - cannot parse arguments! It should be (a)");
            }
            std::string filename_str(filename);

            PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
            Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;

            self->cpp_obj =  new FieldTracker(bx, by,
          	       ax, ay, ex,  epx, l, zi,  zf, ds,  niters, resid,
          	       xrefi, yrefi, eulerai, eulerbi, eulergi, cpp_bunch, filename_str);

            ((FieldTracker*) self->cpp_obj)->setPyWrapper((PyObject*) self);
            return 0;
        }

        /** Performs the collimation tracking of the bunch */
        static PyObject* FieldTracker_trackBunch(PyObject *self, PyObject *args){
            FieldTracker* cpp_FieldTracker = (FieldTracker*)((pyORBIT_Object*) self)->cpp_obj;
            PyObject* pyBunch;
            if(!PyArg_ParseTuple(args,"O:trackBunch",&pyBunch)){
                ORBIT_MPI_Finalize("FieldTracker - trackBunch(Bunch* bunch) - parameter are needed.");
            }
            PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
            if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
                ORBIT_MPI_Finalize("FieldTracker - trackBunch(Bunch* bunch) - method needs a Bunch.");
            }

            Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
            cpp_FieldTracker->trackBunch(cpp_bunch);
            Py_INCREF(Py_None);
            return Py_None;
        }


        static PyObject* FieldTracker_setPathVariable(PyObject *self, PyObject *args){
           int i = 0;
           FieldTracker* cpp_FieldTracker = (FieldTracker*)((pyORBIT_Object*) self)->cpp_obj;
        	if(!PyArg_ParseTuple(args,"i:arguments",  &i))
            {
                error("PyBunch - addParticle - cannot parse arguments! It should be (i)");
            }
        	cpp_FieldTracker->setPathVariable(i);
        	return Py_None;
        }

        //-----------------------------------------------------
        //destructor for python FieldTracker class (__del__ method).
        //-----------------------------------------------------
        static void FieldTracker_del(pyORBIT_Object* self){
            //std::cerr<<"The FieldTracker __del__ has been called!"<<std::endl;
            delete ((FieldTracker*)self->cpp_obj);
            self->ob_base.ob_type->tp_free((PyObject*)self);
        }

        // definition of the methods of the python FieldTrackerwrapper class
        // they will be vailable from python level
        static PyMethodDef FieldTrackerClassMethods[] = {
            { "trackBunch",FieldTracker_trackBunch,METH_VARARGS,"Performs the field tracking of the bunch."},
            {"setPathVariable", FieldTracker_setPathVariable,METH_VARARGS,"Determines whether or not to output path"},
            {NULL}
        };

        static PyMemberDef FieldTrackerClassMembers [] = {
            {NULL}
        };

        //new python FieldTracker wrapper type definition
        static PyTypeObject pyORBIT_FieldTracker_Type = {
            PyVarObject_HEAD_INIT(NULL, 0)
            "FieldTracker", /*tp_name*/
            sizeof(pyORBIT_Object), /*tp_basicsize*/
            0, /*tp_itemsize*/
            (destructor) FieldTracker_del , /*tp_dealloc*/
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
            "The FieldTracker python wrapper", /* tp_doc */
            0, /* tp_traverse */
            0, /* tp_clear */
            0, /* tp_richcompare */
            0, /* tp_weaklistoffset */
            0, /* tp_iter */
            0, /* tp_iternext */
            FieldTrackerClassMethods, /* tp_methods */
            FieldTrackerClassMembers, /* tp_members */
            0, /* tp_getset */
            0, /* tp_base */
            0, /* tp_dict */
            0, /* tp_descr_get */
            0, /* tp_descr_set */
            0, /* tp_dictoffset */
            (initproc) FieldTracker_init, /* tp_init */
            0, /* tp_alloc */
            FieldTracker_new, /* tp_new */
        };


        //--------------------------------------------------
        //Initialization of the pyFieldTracker class
        //--------------------------------------------------
static struct PyModuleDef cModFieldTrack = {
    PyModuleDef_HEAD_INIT,
    "fieldtracker", "FieldTracker class",
    -1,
     FieldTrackerClassMethods
     };
        PyMODINIT_FUNC initfieldtracker(){
            //check that the FieldTracker wrapper is ready
            if (PyType_Ready(&pyORBIT_FieldTracker_Type) < 0) return NULL;
            Py_INCREF(&pyORBIT_FieldTracker_Type);
            //create new module
            PyObject* module = PyModule_Create(&cModFieldTrack);
            PyModule_AddObject(module, "FieldTracker", (PyObject *)&pyORBIT_FieldTracker_Type);
	    return module;
 }

#ifdef __cplusplus
    }
#endif


}
