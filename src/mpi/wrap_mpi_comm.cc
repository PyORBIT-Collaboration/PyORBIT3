#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

namespace wrap_orbit_mpi_comm{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	static PyObject* mpi_comm_new(PyTypeObject *type, PyObject *Py_UNUSED(args), PyObject *Py_UNUSED(kwds))
	{
		pyORBIT_MPI_Comm* self;
		self = (pyORBIT_MPI_Comm *) type->tp_alloc(type, 0);
		self->comm = MPI_COMM_WORLD;
		return (PyObject *) self;
	}

	static int mpi_comm_init(pyORBIT_MPI_Comm *Py_UNUSED(self), PyObject *args, PyObject *Py_UNUSED(kwds)){
		if(PyTuple_Size(args) != 0){
			error("MPI_Comm constructor cannot have an input parameter.");
		}
		return 0;
	}

	static PyObject* mpi_comm_free(PyObject *self, PyObject *Py_UNUSED(ignored)){
		pyORBIT_MPI_Comm* pyMPI_Comm = (pyORBIT_MPI_Comm*) self;
		if(pyMPI_Comm->comm != MPI_COMM_WORLD && pyMPI_Comm->comm != MPI_COMM_SELF){
			ORBIT_MPI_Comm_free(&pyMPI_Comm->comm);
		}
		pyMPI_Comm->comm = MPI_COMM_WORLD;
		Py_INCREF(Py_None);
		return Py_None;
	}

	static void mpi_comm_del(pyORBIT_MPI_Comm* self){
		MPI_Comm comm = self->comm;
		if(comm != MPI_COMM_NULL && comm != MPI_COMM_WORLD && comm != MPI_COMM_SELF){
			ORBIT_MPI_Comm_free(&comm);
		}
		self->ob_base.ob_type->tp_free((PyObject*)self);
	}

	static PyMethodDef MPI_CommClassMethods[] = {
		{ "free",       mpi_comm_free ,METH_NOARGS,"Free MPI communicator."},
		{NULL}
	};

	static PyMemberDef MPI_CommClassMembers[] = {
		{NULL}
	};

	PyTypeObject pyORBIT_MPI_Comm_Type = {
		PyVarObject_HEAD_INIT(NULL, 0)
		"MPI_Comm",
		sizeof(pyORBIT_MPI_Comm),
		0,
		(destructor) mpi_comm_del,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
		"The MPI_Comm python wrapper",
		0,
		0,
		0,
		0,
		0,
		0,
		MPI_CommClassMethods,
		MPI_CommClassMembers,
		0,
		0,
		0,
		0,
		0,
		0,
		(initproc) mpi_comm_init,
		0,
		mpi_comm_new,
	};

	void init_orbit_mpi_comm(PyObject* module){
		if (PyType_Ready(&pyORBIT_MPI_Comm_Type) < 0) return;
		Py_INCREF(&pyORBIT_MPI_Comm_Type);

		PyObject * comm_module = PyModule_New("mpi_comm");
		PyModule_AddObject(comm_module, "MPI_Comm", (PyObject *)&pyORBIT_MPI_Comm_Type);
		Py_INCREF(comm_module);

		pyORBIT_MPI_Comm* pyMPI_Comm_WORLD = PyObject_New(pyORBIT_MPI_Comm,&pyORBIT_MPI_Comm_Type);
		pyMPI_Comm_WORLD->comm = MPI_COMM_WORLD;
		Py_INCREF((PyObject *) pyMPI_Comm_WORLD);

		pyORBIT_MPI_Comm* pyMPI_Comm_SELF = PyObject_New(pyORBIT_MPI_Comm,&pyORBIT_MPI_Comm_Type);
		pyMPI_Comm_SELF->comm = MPI_COMM_SELF;
		Py_INCREF((PyObject *) pyMPI_Comm_SELF);

		pyORBIT_MPI_Comm* pyMPI_Comm_NULL = PyObject_New(pyORBIT_MPI_Comm,&pyORBIT_MPI_Comm_Type);
		pyMPI_Comm_NULL->comm = MPI_COMM_NULL;
		Py_INCREF((PyObject *) pyMPI_Comm_NULL);

		PyModule_AddObject(comm_module, "MPI_COMM_WORLD", (PyObject *) pyMPI_Comm_WORLD);
		PyModule_AddObject(comm_module, "MPI_COMM_SELF", (PyObject *) pyMPI_Comm_SELF);
		PyModule_AddObject(comm_module, "MPI_COMM_NULL", (PyObject *) pyMPI_Comm_NULL);

		PyModule_AddObject(module, "mpi_comm", comm_module);
	}

	pyORBIT_MPI_Comm* newMPI_Comm(){
		pyORBIT_MPI_Comm* pyMPI_Comm = PyObject_New(pyORBIT_MPI_Comm,&pyORBIT_MPI_Comm_Type);
		pyMPI_Comm->comm = MPI_COMM_WORLD;
		Py_INCREF((PyObject *) pyMPI_Comm);
		return pyMPI_Comm;
	}

	void freeMPI_Comm(pyORBIT_MPI_Comm* pyMPI_Comm){
		Py_DECREF(pyMPI_Comm);
	}

#ifdef __cplusplus
}
#endif

}
