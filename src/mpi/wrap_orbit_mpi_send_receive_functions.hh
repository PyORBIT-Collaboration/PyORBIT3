//------------------------------------------------------------------
// Wrappers for MPI functions related to the Send-Receive operations
//------------------------------------------------------------------
static PyObject* mpi_barrier(PyObject *self, PyObject *args){
	PyObject* pyO_comm;
	if(!PyArg_ParseTuple(	args,"O:mpi_barrier",&pyO_comm)){
		error("MPI_ Barrier(MPI_Comm comm)- needs 1 param.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO_comm;
	int res = ORBIT_MPI_Barrier(pyComm->comm);
	if(res != MPI_SUCCESS){
		error("MPI_ Barrier(MPI_Comm comm)- fatal error. STOP.");
	}
	return Py_BuildValue("i",res);
}

static PyObject* mpi_wait(PyObject *self, PyObject *args){
	PyObject* pyO_request; PyObject* pyO_status;
	if(!PyArg_ParseTuple(	args,"OO:mpi_wait",&pyO_request, &pyO_status)){
		error("MPI_ Wait(MPI_Request,MPI_Status)- needs 2 params.");
	}
	pyORBIT_MPI_Request* pyRequest = (pyORBIT_MPI_Request*) pyO_request;
	pyORBIT_MPI_Status* pyStatus = (pyORBIT_MPI_Status*) pyO_status;
	int res = ORBIT_MPI_Wait(&pyRequest->request,&pyStatus->status);
	if(res != MPI_SUCCESS){
		error("MPI_ Wait(MPI_Request,MPI_Status)- fatal error. STOP.");
	}
	return Py_BuildValue("i",res);
}

static PyObject* allreduce_scalar_int(PyObject* obj, pyORBIT_MPI_Op* pyOp, pyORBIT_MPI_Comm* pyComm) {
  long v = PyLong_AsLong(obj);
  if (PyErr_Occurred()) return NULL;

  int val = (int)v;
  ORBIT_MPI_Allreduce(MPI_IN_PLACE, &val, 1, MPI_INT, pyOp->op, pyComm->comm);
  return Py_BuildValue("i", val);
}

static PyObject* allreduce_scalar_double(PyObject* obj, pyORBIT_MPI_Op* pyOp, pyORBIT_MPI_Comm* pyComm) {
  double val = PyFloat_AsDouble(obj);
  if (PyErr_Occurred()) return NULL;

  ORBIT_MPI_Allreduce(MPI_IN_PLACE, &val, 1, MPI_DOUBLE, pyOp->op, pyComm->comm);
  return Py_BuildValue("d", val);
}

static PyObject* allreduce_sequence_int(PyObject* obj, pyORBIT_MPI_Op* pyOp, pyORBIT_MPI_Comm* pyComm) {
  PyObject *seq = PySequence_Fast(obj, "MPI_Allreduce(...) expected a sequence of ints");
  if (seq == NULL) return NULL;

  Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
  PyObject *pyRes = PyTuple_New(size);

  if (pyRes == NULL) {
	Py_DECREF(seq);
	return NULL;
  }

  int buff_index = 0;
  int *buf = BufferStore::getBufferStore()->getFreeIntArr(buff_index, (int)size);

  for (Py_ssize_t i = 0; i < size; ++i) {
    long v = PyLong_AsLong(PySequence_Fast_GET_ITEM(seq, i));
	if (PyErr_Occurred()) {
	  BufferStore::getBufferStore()->setUnusedIntArr(buff_index);
	  Py_DECREF(pyRes);
	  Py_DECREF(seq);
	  return NULL;
	}
	buf[i] = (int)v;
  }

  ORBIT_MPI_Allreduce(MPI_IN_PLACE, buf, (int)size, MPI_INT, pyOp->op, pyComm->comm);

  for (Py_ssize_t i = 0; i < size; ++i) {
	PyObject *item = Py_BuildValue("i", buf[i]);
	if (item == NULL || PyTuple_SetItem(pyRes, i, item) != 0) {
	  Py_XDECREF(item);
	  BufferStore::getBufferStore()->setUnusedIntArr(buff_index);
	  Py_DECREF(pyRes);
	  Py_DECREF(seq);
	  return NULL;
	}
  }

  BufferStore::getBufferStore()->setUnusedIntArr(buff_index);
  Py_DECREF(seq);
  return pyRes;
}

static PyObject* allreduce_sequence_double(PyObject* obj, pyORBIT_MPI_Op* pyOp, pyORBIT_MPI_Comm* pyComm) {
  PyObject *seq = PySequence_Fast(obj, "MPI_Allreduce(...) expected a sequence of ints");
  if (seq == NULL) return NULL;

  Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
  PyObject *pyRes = PyTuple_New(size);

  if (pyRes == NULL) {
	Py_DECREF(seq);
	return NULL;
  }

  int buff_index = 0;
  double *buf = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index, (int)size);

  for (Py_ssize_t i = 0; i < size; ++i) {
	double v = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
	if (PyErr_Occurred()) {
	  BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index);
	  Py_DECREF(pyRes);
	  Py_DECREF(seq);
	  return NULL;
	}
	buf[i] = v;
  }

  ORBIT_MPI_Allreduce(MPI_IN_PLACE, buf, (int)size, MPI_DOUBLE, pyOp->op, pyComm->comm);

  for (Py_ssize_t i = 0; i < size; ++i) {
	PyObject *item = Py_BuildValue("d", buf[i]);
	if (item == NULL || PyTuple_SetItem(pyRes, i, item) != 0) {
	  Py_XDECREF(item);
	  BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index);
	  Py_DECREF(pyRes);
	  Py_DECREF(seq);
	  return NULL;
	}
  }

  BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index);
  Py_DECREF(seq);
  return pyRes;
}

static PyObject* mpi_allreduce(PyObject *self, PyObject *args){
  PyObject *pyO_arr, *pyO_datatype, *pyO_op,  *pyO_comm;

  if(!PyArg_ParseTuple(args,"OOOO:mpi_allreduce", &pyO_arr, &pyO_datatype, &pyO_op, &pyO_comm)) {
  	error("MPI_Allreduce([...data],MPI_Datatype type,MPI_Op op,MPI_Comm out) - needs 4 params.");
  }

  pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO_comm;
  pyORBIT_MPI_Datatype* pyDatatype = (pyORBIT_MPI_Datatype*) pyO_datatype;
  pyORBIT_MPI_Op* pyOp = (pyORBIT_MPI_Op*) pyO_op;

  //check the data type
  MPI_Datatype dtype = pyDatatype->datatype;

  if(dtype != MPI_INT && dtype != MPI_DOUBLE){
  	error("MPI_Allreduce(...)  data type should be INT or DOUBLE. STOP.");
  }

 //check if it is a sequence or a scalar
 if (PySequence_Check(pyO_arr) == 0) {
   if (dtype == MPI_INT) {
     return allreduce_scalar_int(pyO_arr, pyOp, pyComm);
   }
   return allreduce_scalar_double(pyO_arr, pyOp, pyComm);
 } else {
   if (dtype == MPI_INT) {
     return allreduce_sequence_int(pyO_arr, pyOp, pyComm);
   }
   return allreduce_sequence_double(pyO_arr, pyOp, pyComm);
 }
}

static PyObject* mpi_bcast(PyObject *self, PyObject *args){
	PyObject* pyO_arr; PyObject* pyO_datatype; PyObject* pyO_comm;
	int rank;
	if(!PyArg_ParseTuple(	args,"OOiO:mpi_bcast",&pyO_arr,&pyO_datatype,&rank,&pyO_comm)){
		error("MPI_Bcast([...data],MPI_Datatype type,int rank,MPI_Comm out) - needs 4 params.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO_comm;
	pyORBIT_MPI_Datatype* pyDatatype = (pyORBIT_MPI_Datatype*) pyO_datatype;
	int rank_local = -1;
	ORBIT_MPI_Comm_rank(pyComm->comm,&rank_local);
	//check if it is not sequence
	int is_seq = 0;
	if(rank_local == rank && PySequence_Check(pyO_arr) == 1){
		is_seq = 1;
	}
	ORBIT_MPI_Bcast(&is_seq,1,MPI_INT,rank,pyComm->comm);
	//it is NOT SEQUENCE
	if(is_seq == 0){
		if(pyDatatype->datatype == MPI_INT){
			int val = 0;
			if(rank_local == rank){
				val = (int) PyLong_AsLong(pyO_arr);
			}
			ORBIT_MPI_Bcast(&val,1,MPI_INT,rank,pyComm->comm);
			return Py_BuildValue("i",val);
		}
		if(pyDatatype->datatype == MPI_DOUBLE){
			double val = 0;
			if(rank_local == rank){
				val = PyFloat_AsDouble(pyO_arr);
			}
			ORBIT_MPI_Bcast(&val,1,MPI_DOUBLE,rank,pyComm->comm);
			return Py_BuildValue("d",val);
		}
		error("MPI_Bcast(...) - use only INT and DOUBLE data type as scalar");
	}
	//it IS A SEQUENCE
  int size = 0;
	if(rank_local == rank){
		size =  PySequence_Size(pyO_arr);
	}
	ORBIT_MPI_Bcast(&size,1,MPI_INT,rank,pyComm->comm);
	PyObject* pyRes = Py_None;
	//data is an INT array
	if(pyDatatype->datatype == MPI_INT){
		pyRes = PyTuple_New(size);
		int buff_index = 0;
		int* arr =  BufferStore::getBufferStore()->getFreeIntArr(buff_index,size);
		if(rank_local == rank){
			for(int i = 0; i < size; i++){
				arr[i]= (int) PyLong_AsLong(PySequence_Fast_GET_ITEM(pyO_arr, i));
				if(arr[i] < 0){
					error("MPI_Bcast([...data],MPI_Datatype type,int rank,MPI_Comm out) [...data] is not good.");
				}
			}
		}
		ORBIT_MPI_Bcast(arr,size,MPI_INT,rank,pyComm->comm);
		for(int i = 0; i < size; i++){
			if(PyTuple_SetItem(pyRes,i,Py_BuildValue("i",arr[i])) != 0){
				error("MPI_Bcast(...)  cannot create a resulting tuple.");
			}
		}
		BufferStore::getBufferStore()->setUnusedIntArr(buff_index);
	}
	//data is an DOUBLE array
	if(pyDatatype->datatype == MPI_DOUBLE){
		pyRes = PyTuple_New(size);
		int buff_index = 0;
		double* arr =  BufferStore::getBufferStore()->getFreeDoubleArr(buff_index,size);
		if(rank_local == rank){
			for(int i = 0; i < size; i++){
				arr[i]= PyFloat_AsDouble(PySequence_Fast_GET_ITEM(pyO_arr, i));
			}
		}
		ORBIT_MPI_Bcast(arr,size,MPI_DOUBLE,rank,pyComm->comm);
		for(int i = 0; i < size; i++){
			if(PyTuple_SetItem(pyRes,i,Py_BuildValue("d",arr[i])) != 0){
				error("MPI_Bcast(...)  cannot create a resulting tuple.");
			}
		}
		BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index);
	}
	//data is a pyString
	if(pyDatatype->datatype == MPI_CHAR){
		char* arr = NULL;
		int buff_index = -1;
		if(rank_local == rank){
			arr = (char *)PyUnicode_AsUTF8(pyO_arr);
		} else {
			arr = BufferStore::getBufferStore()->getFreeCharArr(buff_index,size);
		}
		if(arr == NULL){
			error("MPI_Bcast(...)  data type could be INT, DOUBLE, or CHAR(string). STOP.");
		}
		ORBIT_MPI_Bcast(arr,size,MPI_CHAR,rank,pyComm->comm);
		pyRes = Py_BuildValue("s#",arr,size);
		if(buff_index >= 0){
			BufferStore::getBufferStore()->setUnusedCharArr(buff_index);
		}
	}
	if(pyRes == Py_None){
		error("MPI_Bcast(...)  data type could be INT, DOUBLE, or CHAR(string). STOP.");
	}

	return pyRes;
}

static PyObject* mpi_send(PyObject *self, PyObject *args){
	PyObject* pyO_arr; PyObject* pyO_datatype; PyObject* pyO_comm;
	int dest, tag;
	if(!PyArg_ParseTuple(	args,"OOiiO:mpi_send",&pyO_arr,&pyO_datatype,&dest,&tag,&pyO_comm)){
		error("MPI_Send([...data],MPI_Datatype type,int dest, int tag, MPI_Comm) - needs 5 params.");
	}
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO_comm;
	pyORBIT_MPI_Datatype* pyDatatype = (pyORBIT_MPI_Datatype*) pyO_datatype;
	int rank = -1;
	ORBIT_MPI_Comm_rank(pyComm->comm, &rank);
	//check if it is not sequence
	if(PySequence_Check(pyO_arr) != 1){
		if(pyDatatype->datatype == MPI_INT){
			int val = (int) PyLong_AsLong(pyO_arr);
			int res = ORBIT_MPI_Send(&val,1,MPI_INT,dest,tag,pyComm->comm);
			if(res != MPI_SUCCESS){
				std::cerr << "Error MPI_Send(...) rank="<< rank <<" dest_rank="<<dest<<" tag="<<tag<<std::endl;
				error("MPI_Send(...) - use fatal error. Stop");
			}
			return Py_BuildValue("i",res);
		}
		if(pyDatatype->datatype == MPI_DOUBLE){
			double val = PyFloat_AsDouble(pyO_arr);
			int res = ORBIT_MPI_Send(&val,1,MPI_DOUBLE,dest,tag,pyComm->comm);
			if(res != MPI_SUCCESS){
				std::cerr << "Error MPI_Send(...) rank="<< rank <<" dest_rank="<<dest<<" tag="<<tag<<std::endl;
				error("MPI_Send(...) - use fatal error. Stop");
			}
			return Py_BuildValue("i",res);
		}
		error("MPI_Send(...) - use only INT and DOUBLE data type as scalar");
	}
	//it IS A SEQUENCE
	int res = 0;
  int size = PySequence_Size(pyO_arr);
	int buff_index = -1;
	//data is an INT array
	if(pyDatatype->datatype == MPI_INT){
		int* arr =  BufferStore::getBufferStore()->getFreeIntArr(buff_index,size);
		for(int i = 0; i < size; i++){
			arr[i]= (int) PyLong_AsLong(PySequence_Fast_GET_ITEM(pyO_arr, i));
		}
		res = ORBIT_MPI_Send(arr,size,MPI_INT,dest,tag,pyComm->comm);
		BufferStore::getBufferStore()->setUnusedIntArr(buff_index);
	}
	//data is an DOUBLE array
	if(pyDatatype->datatype == MPI_DOUBLE){
		double* arr =  BufferStore::getBufferStore()->getFreeDoubleArr(buff_index,size);
		for(int i = 0; i < size; i++){
			arr[i]= PyFloat_AsDouble(PySequence_Fast_GET_ITEM(pyO_arr, i));
		}
		res = ORBIT_MPI_Send(arr,size,MPI_DOUBLE,dest,tag,pyComm->comm);
		BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index);
	}
	//data is a pyString
	if(pyDatatype->datatype == MPI_CHAR){
		char* arr = (char *)PyUnicode_AsUTF8(pyO_arr);
		res = ORBIT_MPI_Send(arr,size,MPI_CHAR,dest,tag,pyComm->comm);
	}
	if(res != MPI_SUCCESS){
		error("MPI_Send(...)  data type could be INT[], DOUBLE[], or CHAR(string). STOP.");
	}
	return Py_BuildValue("i",res);
}

static PyObject* mpi_recv(PyObject *self, PyObject *args){
	PyObject* pyO_datatype; PyObject* pyO_comm;
	int source, tag;
	if(!PyArg_ParseTuple(	args,"OiiO:mpi_recv",&pyO_datatype,&source,&tag,&pyO_comm)){
		error("MPI_Recv(MPI_Datatype type,int source, int tag,MPI_Comm comm) - needs 4 params.");
	}
	pyORBIT_MPI_Datatype* pyDatatype = (pyORBIT_MPI_Datatype*) pyO_datatype;
	pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) pyO_comm;
	MPI_Status status;
	int rank = -1;
	ORBIT_MPI_Comm_rank(pyComm->comm,&rank);
	//check the size of the message
	ORBIT_MPI_Probe(source,tag,pyComm->comm,&status);
	int size;
	ORBIT_MPI_Get_count(&status, pyDatatype->datatype, &size);
	int res = 0;
	if(size < 1){
		std::cerr << "Error MPI_Recv(...) rank="<< rank <<" source_rank="<<source<<" tag="<<tag<<std::endl;
		error("MPI_Recv(...) - fatal error. STOP.");
	}
	if(size == 1){
		if(pyDatatype->datatype == MPI_INT){
			int val = 0;
			res = ORBIT_MPI_Recv(&val,1,MPI_INT,source,tag,pyComm->comm,&status);
			if(res != MPI_SUCCESS){
				error("MPI_Recv(...) - fatal error. STOP.");
			}
			return Py_BuildValue("i",val);
		}
		if(pyDatatype->datatype == MPI_DOUBLE){
		  double val = 0;
			res = ORBIT_MPI_Recv(&val,1,MPI_DOUBLE,source,tag,pyComm->comm,&status);
			if(res != MPI_SUCCESS){
				error("MPI_Recv(...) - fatal error. STOP.");
			}
			return Py_BuildValue("d",val);
		}
		if(pyDatatype->datatype == MPI_CHAR){
		  char val;
			res = ORBIT_MPI_Recv(&val,1,MPI_CHAR,source,tag,pyComm->comm,&status);
			if(res != MPI_SUCCESS){
				error("MPI_Recv(...) - fatal error. STOP.");
			}
			return Py_BuildValue("s#",&val,1);
		}
		error("MPI_Recv(...) - data could be INT, DOUBLE, or CHAR(string). STOP.");
	}
	PyObject* pyRes =  Py_None;
	if(pyDatatype->datatype == MPI_INT){
		pyRes = PyTuple_New(size);
		int buff_index = 0;
		int* arr =  BufferStore::getBufferStore()->getFreeIntArr(buff_index,size);
		res = ORBIT_MPI_Recv(arr,size,MPI_INT,source,tag,pyComm->comm,&status);
		if(res != MPI_SUCCESS){
			error("MPI_Recv(...) - fatal error. STOP.");
		}
		for(int i = 0; i < size; i++){
			if(PyTuple_SetItem(pyRes,i,Py_BuildValue("i",arr[i])) != 0){
				error("MPI_Recv(...)  cannot create a resulting tuple.");
			}
		}
		BufferStore::getBufferStore()->setUnusedIntArr(buff_index);
	}
	if(pyDatatype->datatype == MPI_DOUBLE){
		pyRes = PyTuple_New(size);
		int buff_index = 0;
		double* arr =  BufferStore::getBufferStore()->getFreeDoubleArr(buff_index,size);
		res = ORBIT_MPI_Recv(arr,size,MPI_DOUBLE,source,tag,pyComm->comm,&status);
		if(res != MPI_SUCCESS){
			error("MPI_Recv(...) - fatal error. STOP.");
		}
		for(int i = 0; i < size; i++){
			if(PyTuple_SetItem(pyRes,i,Py_BuildValue("d",arr[i])) != 0){
				error("MPI_Recv(...)  cannot create a resulting tuple.");
			}
		}
		BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index);
	}
	if(pyDatatype->datatype == MPI_CHAR){
		int buff_index = 0;
		char* arr =  BufferStore::getBufferStore()->getFreeCharArr(buff_index,size);
		res = ORBIT_MPI_Recv(arr,size,MPI_CHAR,source,tag,pyComm->comm,&status);
		if(res != MPI_SUCCESS){
			error("MPI_Recv(...) - fatal error. STOP.");
		}
		pyRes = Py_BuildValue("s#",arr,size);
		BufferStore::getBufferStore()->setUnusedCharArr(buff_index);
	}
	if(pyRes == Py_None){
		error("MPI_Recv(...)  data type could be INT, DOUBLE, or CHAR(string). STOP.");
	}
	return pyRes;
}
