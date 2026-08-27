//////////////////////////////// -*- C++ -*- //////////////////////////////
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "wrap_bunch.hh"
// #include "modsupport.h"
// #include "pyerrors.h"
#include "wrap_syncpart.hh"
#include "wrap_bunch_twiss_analysis.hh"
#include "wrap_bunch_tune_analysis.hh"
#include "wrap_synch_part_redefinition_z_de.hh"

#include "pyORBIT_Object.hh"

#ifdef PyORBIT_EXPERIMENTAL_WITH_NUMPY
#include <numpy/arrayobject.h>

static int ensure_numpy() {
  static int numpy_initialized = 0;
  if (!numpy_initialized) {
    import_array1(-1);
    numpy_initialized = 1;
  }
  return 0;
}
#endif // PyORBIT_EXPERIMENTAL_WITH_NUMPY

#include "Bunch.hh"
#include "ParticleAttributesFactory.hh"

namespace wrap_orbit_bunch{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }
    //---------------------------------------------------------
    //Python Bunch class definition
    //---------------------------------------------------------

    //constructor for python class wrapping Bunch instance
    //It never will be called directly
    static PyObject* Bunch_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        return (PyObject *) self;
    }

  //initializator for python Bunch class
  //this is implementation of the __init__ method
  static int Bunch_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		//std::cerr<<"The Bunch __init__ has been called!"<<std::endl;
		//instantiation of a new c++ Bunch
		self->cpp_obj = (void*) new Bunch();
		((Bunch*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//This is the way to create new class instance from the C-level
		// Template: PyObject* PyObject_CallMethod(	PyObject *o, char *method, char *format, ...)
		//see Python/C API documentation
		//It will create a SyncParticle object and set the reference to it from pyBunch
		PyObject* mod = PyImport_ImportModule("orbit.core.bunch");
		PyObject* pySyncPart = PyObject_CallMethod(mod,const_cast<char*>("SyncParticle"),const_cast<char*>("O"),self);

		//the references should be decreased because they were created as "new reference"
		Py_DECREF(pySyncPart);
		Py_DECREF(mod);
    return 0;
  }

  //---------------------------------------------------------------
  //
  // methods related to synchronous particle etc.
  //
  //----------------------------------------------------------------

  //returns the SyncPart python class wrapper instance
    static PyObject* Bunch_getSyncParticle(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
        PyObject* pySyncPart = cpp_bunch->getSyncPart()->getPyWrapper();
        Py_INCREF(pySyncPart);
    return pySyncPart;
  }

  //---------------------------------------------------------------
  //
  // set and get MPI Communicators
  //
  //----------------------------------------------------------------

    //returns the local MPI Comm for this bunch
    static PyObject* Bunch_getMPIComm(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
        PyObject* pyMPIComm = (PyObject*) cpp_bunch->getMPI_Comm_Local();
        Py_INCREF(pyMPIComm);
    return pyMPIComm;
  }

    //sets a new local MPI Comm for this bunch
    static PyObject* Bunch_setMPIComm(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
        cpp_bunch->setMPI_Comm_Local( (pyORBIT_MPI_Comm*) arg);
        Py_INCREF(Py_None);
    return Py_None;
  }

  //---------------------------------------------------------------
  //
  // add and remove particles, compress etc.
  //
  //----------------------------------------------------------------

  //adds a particle to the Bunch object
  //this is implementation of the  addParticle(...) method
  static PyObject* Bunch_addParticle(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;

    double x = 0.;  double xp = 0.; double y = 0.;
    double yp = 0.; double z = 0.;  double zp = 0.;

    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(args,"dddddd:coordinates",&x,&xp,&y,&yp,&z,&zp)){
      error("PyBunch - addParticle - cannot parse arguments! It should be (x,xp,y,yp,z,zp)");
    }
    int ind = cpp_bunch->addParticle(x,xp,y,yp,z,zp);
    return Py_BuildValue("i",ind);
  }


  //removes a particle to the Bunch object
  //returns the number of particles in the bunch
  //this is implementation of the deleteParticle(int index)  method
  static PyObject* Bunch_deleteParticle(PyObject *self, PyObject *arg){
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    int ind;

    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"i:deleteParticle",&ind)){
      return NULL;
    }

    cpp_bunch->deleteParticle(ind);
    int size = cpp_bunch->getSize();

    return Py_BuildValue("i",size);
  }

  //removes a particle to the Bunch object
  //returns the index of removed macro-particle
  //this is implementation of the deleteParticleFast(int index)  method
  static PyObject* Bunch_deleteParticleFast(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    int ind;

    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"i:deleteParticleFast",&ind)){
      return NULL;
    }

    cpp_bunch->deleteParticleFast(ind);
    return Py_BuildValue("i",ind);
  }

  static PyObject* Bunch_recoverParticle(PyObject *self, PyObject *arg){
    Bunch *cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    int ind;

    if(!PyArg_Parse(arg,"i:recoverParticle",&ind)){
      return NULL;
    }

    cpp_bunch->recoverParticle(ind);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //removes all particles from the Bunch object
  //this is implementation of the deleteAllParticles()  method
  static PyObject* Bunch_deleteAllParticles(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->deleteAllParticles();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //compress the bunch. This method should be called after deleting one
  //  or more macro-particles
  //this is implementation of the compress()  method
  static PyObject* Bunch_compress(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->compress();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //---------------------------------------------------------------
  //
  // related to the macro-particles' coordinates
  //
  //----------------------------------------------------------------

  //Sets or returns x coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns x-coordinate
  //  (index, value) - sets the new value to the x-coordinate
  //this is implementation of the x(int index)  method
  static PyObject* Bunch_x(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"i:x",&index)){
          error("PyBunch - x(index) - index is needed");
        }
        val = cpp_bunch->x(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"id:x",&index,&val)){
          error("PyBunch - x(index, value) - index and value are needed");
        }
        cpp_bunch->x(index) = val;
      }
            return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call x(index) or x(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns y coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns y-coordinate
  //  (index, value) - sets the new value to the y-coordinate
  //this is implementation of the y(int index)  method
  static PyObject* Bunch_y(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"i:y",&index)){
          error("PyBunch - y(index) - index is needed");
        }
        val = cpp_bunch->y(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"id:y",&index,&val)){
          error("PyBunch - y(index, value) - index and value are needed");
        }
        cpp_bunch->y(index) = val;
      }
            return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call y(index) or y(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns z coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns z(phi)-coordinate
  //  (index, value) - sets the new value to the z(phi)-coordinate
  //this is implementation of the (z or phi)(int index)  method
  static PyObject* Bunch_z(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"i:z",&index)){
          error("PyBunch - z(index) - index is needed");
        }
        val = cpp_bunch->z(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"id:z",&index,&val)){
          error("PyBunch - z(index, value) - index and value are needed");
        }
        cpp_bunch->z(index) = val;
      }
            return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call z(index) or z(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }


  //Sets or returns px coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns px-coordinate
  //  (index, value) - sets the new value to the px-coordinate
  //this is implementation of the px(int index)  method
  static PyObject* Bunch_px(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"i:px",&index)){
          error("PyBunch - px(index) - index is needed");
        }
        val = cpp_bunch->px(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"id:px",&index,&val)){
          error("PyBunch - px(index, value) - index and value are needed");
        }
        cpp_bunch->px(index) = val;
      }
            return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call px(index) or px(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns py coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns y-coordinate
  //  (index, value) - sets the new value to the py-coordinate
  //this is implementation of the py(int index)  method
  static PyObject* Bunch_py(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"i:py",&index)){
          error("PyBunch - py(index) - index is needed");
        }
        val = cpp_bunch->py(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"id:py",&index,&val)){
          error("PyBunch - py(index, value) - index and value are needed");
        }
        cpp_bunch->py(index) = val;
      }
            return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call py(index) or py(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns pz or dE coordinate of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns pz(dE)-coordinate
  //  (index, value) - sets the new value to the pz(dE)-coordinate
  //this is implementation of the (pz or dE)(int index)  method
  static PyObject* Bunch_pz(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 get coordinate
    //if nVars == 2 set coordinate
    int nVars = PyTuple_Size(args);

    int index = 0;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"i:pz",&index)){
          error("PyBunch - pz(index) - index is needed");
        }
        val = cpp_bunch->pz(index);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"id:pz",&index,&val)){
          error("PyBunch - pz(index, value) - index and value are needed");
        }
        cpp_bunch->pz(index) = val;
      }
            return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call pz(index) or pz(index,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns flag of the macro-particle
  //  the action is depended on the number of arguments
  //  (index) - returns flag
  //this is implementation of the flag(int index)  method
  static PyObject* Bunch_flag(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    int index = 0;
    if(!PyArg_Parse(arg,"i:flag",&index)){
      return NULL;
    }
    int flag = cpp_bunch->flag(index);
    return Py_BuildValue("i",flag);
  }

    //Wraps long. coords in the bunch
    //ringwrap(ring_length)
  static PyObject* Bunch_ringwrap(PyObject *self, PyObject *arg) {
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    double ring_length = 0.;

        //NO NEW OBJECT CREATED BY PyArg_Parse! //NO NEED OF Py_DECREF()
        if(!PyArg_Parse(arg,"d:ringwrap",&ring_length)){
            return NULL;
        }

        cpp_bunch->ringwrap(ring_length);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //---------------------------------------------------------------
  //
  // related to the bunch predefined attributes
  //
  //----------------------------------------------------------------

  //Sets or returns mass of the macro-particle in MeV
  //  the action is depended on the number of arguments
  //  mass() - returns mass
  //  mass(value) - sets the new value
  //this is implementation of the getMass() and setMass  methods of the Bunch class
  static PyObject* Bunch_mass(PyObject *self, PyObject *args){
        Bunch* bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;

        if(0 == PyTuple_GET_SIZE(args)) { 
            return PyFloat_FromDouble(bunch->getMass());
        }

        double value;
        if (!PyArg_ParseTuple(args, "d:mass", &value)) {
            return NULL;
        }

        bunch->setMass(value);
        return PyFloat_FromDouble(value);
  }

  //Returns classicalRadius of the particle in meters
  static PyObject* Bunch_classicalRadius(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
        double val = cpp_bunch->getClassicalRadius();
        return Py_BuildValue("d",val);
  }

  //Returns B_Rho of the particle in [Tesla*meter]. Parameter is used in TEAPOT
  static PyObject* Bunch_B_Rho(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
        double val = cpp_bunch->getB_Rho();
        return Py_BuildValue("d",val);
  }

  //Sets or returns charge of the macro-particle in e-charge
  //  the action is depended on the number of arguments
  //  charge() - returns charge
  //  charge(value) - sets the new value
  //this is implementation of the getCharge() and setCharge  methods of the Bunch class
  static PyObject* Bunch_charge(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 0 get charge
    //if nVars == 1 set charge
    int nVars = PyTuple_Size(args);

    double val = 0.;

    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        val = cpp_bunch->getCharge();
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"d:charge",&val)){
          error("PyBunch - charge(value) - value is needed");
        }
        cpp_bunch->setCharge(val);
      }
            return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call charge() or charge(value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
    }

  //Sets or returns macroSize of the macro-particle
  //  the action is depended on the number of arguments
  //  macroSize() - returns macroSize
  //  macroSize(value) - sets the new value
  //this is implementation of the getMacroSize() and setMacroSize  methods of the Bunch class
  static PyObject* Bunch_macroSize(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 0 get macroSize
    //if nVars == 1 set macroSize
    int nVars = PyTuple_Size(args);

    double val = 0.;

    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        val = cpp_bunch->getMacroSize();
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"d:macroSize",&val)){
          error("PyBunch - macroSize(value) - value is needed");
        }
        cpp_bunch->setMacroSize(val);
      }
            return Py_BuildValue("d",val);
    }
    else{
      error("PyBunch. You should call macroSize() or macroSize(value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //---------------------------------------------------------------
  //
  // related to the bunch attributes
  //
  //----------------------------------------------------------------

  //initilizes bunch attributes from the bunch file
  static PyObject* Bunch_initBunchAttr(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* file_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"s:initBunchAttr",&file_name)){
      return NULL;
    }
    cpp_bunch->initBunchAttributes(file_name);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns a double bunch attribute
  //  the action is depended on the number of arguments
  //  (attr_name) - returns double-value
  //  (index, value) - sets the new value to the attribute
  //this is implementation of
  // getBunchAttributeDouble(name)
  // setBunchAttributeDouble(name,value) Bunch methods
  static PyObject* Bunch_bunchAttrDouble(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 this is get attribute
    //if nVars == 2 this is set attribute
    int nVars = PyTuple_Size(args);

    const char* attr_name = NULL;
    double val = 0.;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"s:bunchAttrDouble",&attr_name)){
          error("PyBunch - bunchAttrDouble(name) - name are needed");
        }
        std::string attr_name_str(attr_name);
        val = cpp_bunch->getBunchAttributeDouble(attr_name_str);
      return Py_BuildValue("d",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"sd:bunchAttrDouble",&attr_name,&val)){
          error("PyBunch - bunchAttrDouble(name,value) - name and double value are needed");
        }

        std::string attr_name_str(attr_name);
        cpp_bunch->setBunchAttribute( attr_name_str, val);
                return Py_BuildValue("d",val);
      }
    }
    else{
      error("PyBunch. You should call bunchAttrDouble(name) or bunchAttrDouble(name,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Sets or returns a integer bunch attribute
  //  the action is depended on the number of arguments
  //  (attr_name) - returns int-value
  //  (index, value) - sets the new value to the attribute
  //this is implementation of
  // getBunchAttributeInt(name)
  // setBunchAttributeInt(name,value) Bunch methods
  static PyObject* Bunch_bunchAttrInt(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 this is get attribute
    //if nVars == 2 this is set attribute
    int nVars = PyTuple_Size(args);

    const char* attr_name = NULL;
    int val = 0;

    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple!NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"s:bunchAttrInt",&attr_name)){
          error("PyBunch - bunchAttrInt(name) - pyBunch object and name are needed");
        }
        std::string attr_name_str(attr_name);
        val = cpp_bunch->getBunchAttributeInt(attr_name_str);
        return Py_BuildValue("i",val);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"si:bunchAttrInt",&attr_name,&val)){
          error("PyBunch - bunchAttrInt(name,value) - name, and double value are needed");
        }

        std::string attr_name_str(attr_name);
        cpp_bunch->setBunchAttribute( attr_name_str, val);
                return Py_BuildValue("i",val);
      }
    }
    else{
      error("PyBunch. You should call bunchAttrInt(name) or bunchAttrInt(name,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns a list (tuple) of  ther double bunch attribute names
  static PyObject* Bunch_bunchAttrDoubleNames(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    std::vector<std::string> names;
    cpp_bunch->getDoubleBunchAttributeNames(names);
        //create tuple with names
        PyObject* resTuple = PyTuple_New(names.size());
        for(int i = 0, n = names.size(); i < n; i++){
            PyObject* py_nm = PyUnicode_FromString(names[i].c_str());
            if(PyTuple_SetItem(resTuple,i,py_nm)){
                error("PyBunch - bunchAttrDoubleNames - cannot create tuple with bunch attr names");
            }
        }
    return resTuple;
  }

  //Returns a list (tuple) of  ther integer bunch attribute names
  static PyObject* Bunch_bunchAttrIntNames(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    std::vector<std::string> names;
    cpp_bunch->getIntBunchAttributeNames(names);
        //create tuple with names
        PyObject* resTuple = PyTuple_New(names.size());
        for(int i = 0, n = names.size(); i < n; i++){
            PyObject* py_nm = PyUnicode_FromString(names[i].c_str());
            if(PyTuple_SetItem(resTuple,i,py_nm)){
                error("PyBunch - bunchAttrIntNames() - cannot create tuple with bunch attr names");
            }
        }
    return resTuple;
  }

  //Returns 0 or 1. The result is 1 if the bunch has an attribute with a particular name
  static PyObject* Bunch_hasBunchAttrDouble(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"s:hasBunchAttrDouble",&attr_name)){
      return NULL;
    }
    std::string attr_name_str(attr_name);
    int res = cpp_bunch->getBunchAttributes()->hasDoubleAttribute(attr_name_str);
    return Py_BuildValue("i",res);
  }

  //Returns 0 or 1. The result is 1 if the bunch has an attribute with a particular name
  static PyObject* Bunch_hasBunchAttrInt(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"s:hasBunchAttrInt",&attr_name)){
      return NULL;
    }
    std::string attr_name_str(attr_name);
    int res = cpp_bunch->getBunchAttributes()->hasIntAttribute(attr_name_str);
    return Py_BuildValue("i",res);
  }

  //---------------------------------------------------------------
  //
  // related to particles' attributes
  //
  //----------------------------------------------------------------

  //Adds a particles' attributes with a particular name to the bunch
  static PyObject* Bunch_addPartAttr(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
        PyObject* py_attrParamsDict = NULL;
    //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
    if(!PyArg_ParseTuple(args,"s|O:addPartAttr",&attr_name,&py_attrParamsDict)){
      error("PyBunch - addPartAttr(name, [param_dict]) - a particle attr. name are needed");
    }
    std::string attr_name_str(attr_name);
        std::map<std::string,double> part_attr_dict;
        if(py_attrParamsDict != NULL){
            if(!PyDict_Check(py_attrParamsDict)){
                error("PyBunch - addPartAttr(name, [param_dict]) - param_dict is not a dictionary");
            }
            PyObject *key, *value;
            Py_ssize_t pos = 0;
            while (PyDict_Next(py_attrParamsDict, &pos, &key, &value)) {
                if(!PyUnicode_Check(key) || !PyNumber_Check(value)){
                    error("PyBunch - addPartAttr(name, [param_dict]) - param_dict is not a {str:val} dictionary");
                }
                std::string par_name((char *)PyUnicode_AsUTF8(key));
                double d_val = PyFloat_AsDouble(value);
                part_attr_dict[par_name] = d_val;
            }
        }
    cpp_bunch->addParticleAttributes(attr_name_str,part_attr_dict);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Removes a particles' attributes with a particular name from the bunch
  static PyObject* Bunch_removePartAttr(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"s:removePartAttr",&attr_name)){
      return NULL;
    }
    std::string attr_name_str(attr_name);
    cpp_bunch->removeParticleAttributes(attr_name_str);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Removes all particles' attributes from the bunch
  static PyObject* Bunch_removeAllPartAttr(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->removeAllParticleAttributes();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns a list (tuple) of  the particles' attributes names
  static PyObject* Bunch_getPartAttrNames(PyObject *self, PyObject *Py_UNUSED(ignored)){
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    std::vector<std::string> names;
    cpp_bunch->getParticleAttributesNames(names);
   //create tuple with names
   PyObject* resTuple = PyTuple_New(names.size());
   for(int i = 0, n = names.size(); i < n; i++){
     PyObject* py_nm = PyUnicode_FromString(names[i].c_str());
     if(PyTuple_SetItem(resTuple,i,py_nm)){
       error("PyBunch - getPartAttrNames() - cannot create tuple with bunch attr names");
     }
   }
    return resTuple;
  }

    //Returns a dict{"part. attribute name":dict{"key":val}}
  static PyObject* Bunch_getPartAttrDicts(PyObject *self, PyObject *Py_UNUSED(ignored)){
    Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    std::vector<std::string> names;
    cpp_bunch->getParticleAttributesNames(names);
        //make dict{"part. attribute name":dict{"key":val}}
        PyObject* resDict = PyDict_New();
        for(int i = 0, n = names.size(); i < n; i++){
            PyObject* py_param_dict = PyDict_New();
            PyDict_SetItemString(resDict,names[i].c_str(),py_param_dict);
            std::map<std::string,double> param_dict = cpp_bunch->getParticleAttributes(names[i])->parameterDict;
            std::map<std::string,double>::iterator pos;
            for (pos = param_dict.begin(); pos != param_dict.end(); ++pos) {
                    std::string key = pos->first;
                    double val = pos->second;
                    PyDict_SetItemString(py_param_dict,key.c_str(),Py_BuildValue("d",val));
            }
        }
    return resDict;
  }

  //Returns a list (tuple) of the possible particles' attributes names
     static PyObject* Bunch_getPossiblePartAttrNames(PyObject *self, PyObject *Py_UNUSED(ignored)){
         std::vector<std::string> names;
         ParticleAttributesFactory::getParticleAttributesNames(names);
         //create tuple with names
         PyObject* resTuple = PyTuple_New(names.size());
         for(int i = 0, n = names.size(); i < n; i++){
             PyObject* py_nm = PyUnicode_FromString(names[i].c_str());
             if(PyTuple_SetItem(resTuple,i,py_nm)){
                 error("PyBunch - getPossiblePartAttrNames - cannot create tuple with bunch attr names");
             }
         }
         return resTuple;
     }

  //temporary removes and memorizes all particles' attributes names
  static PyObject* Bunch_clearAllPartAttrAndMemorize(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->clearAllParticleAttributesAndMemorize();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //restores all particles' attributes names from memory
  static PyObject* Bunch_restoreAllPartAttrFromMemory(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    cpp_bunch->restoreAllParticleAttributesFromMemory();
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns 0 or 1. The result is 1 if the bunch has a particles' attributes with a particular name
  static PyObject* Bunch_hasPartAttr(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"s:hasPartAttr",&attr_name)){
      return NULL;
    }
    std::string attr_name_str(attr_name);
    int res = cpp_bunch->hasParticleAttributes(attr_name_str);
    return Py_BuildValue("i",res);
  }

  //Returns a list (tuple) of  their bunch particles attribute names specified in the bunch file
  static PyObject* Bunch_readPartAttrNames(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* file_name = NULL;
    std::vector<std::string> names;
        std::map<std::string,std::map<std::string,double> > part_attr_dicts;
    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"s:readPartAttrNames",&file_name)){
      return NULL;
    }
    cpp_bunch->readParticleAttributesNames(file_name,names,part_attr_dicts);
        //create tuple with names
        PyObject* resTuple = PyTuple_New(names.size());
        for(int i = 0, n = names.size(); i < n; i++){
            PyObject* py_nm = PyUnicode_FromString(names[i].c_str());
            if(PyTuple_SetItem(resTuple,i,py_nm)){
                error("PyBunch - readPartAttrNames(fileName) - cannot create tuple with particles attr. names");
            }
        }
    return resTuple;
  }

  //Returns a dictionary with the bunch particles attribute names as keys and
    //dictionaries with parameter:value for each attribute
  static PyObject* Bunch_readPartAttrDicts(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* file_name = NULL;
    std::vector<std::string> names;
        std::map<std::string,std::map<std::string,double> > part_attr_dicts;
    if(!PyArg_Parse(arg,"s:readPartAttrDicts",&file_name)){
      return NULL;
    }
    cpp_bunch->readParticleAttributesNames(file_name,names,part_attr_dicts);
        PyObject* resDict = PyDict_New();
        for(int i = 0, n = names.size(); i < n; i++){
            if(part_attr_dicts.count(names[i]) > 0){
                PyObject* py_param_dict = PyDict_New();
                PyDict_SetItemString(resDict,names[i].c_str(),py_param_dict);
                std::map<std::string,double> param_dict = part_attr_dicts[names[i]];
                std::map<std::string,double>::iterator pos;
                for (pos = param_dict.begin(); pos != param_dict.end(); ++pos) {
                    std::string key = pos->first;
                    double val = pos->second;
                    PyDict_SetItemString(py_param_dict,key.c_str(),Py_BuildValue("d",val));
                }
            }
        }
        return resDict;
  }

  //initilizes particles' attributes from the bunch file
  static PyObject* Bunch_readPartAttr(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* file_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"s:readPartAttr",&file_name)){
      return NULL;
    }
    cpp_bunch->readParticleAttributes(file_name);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Returns the number of variables in the particles' attributes with a particular name
  static PyObject* Bunch_getPartAttrSize(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    const char* attr_name = NULL;
    //NO NEW OBJECT CREATED BY PyArg_Parse! NO NEED OF Py_DECREF()
    if(!PyArg_Parse(arg,"s:getPartAttrSize",&attr_name)){
      return NULL;
    }
    std::string attr_name_str(attr_name);
    int size = cpp_bunch->getParticleAttributes(attr_name_str)->getAttSize();
    return Py_BuildValue("i",size);
  }

  //Sets or returns a particles' attributes' value
  //  the action is depended on the number of arguments
  //  (attr_name,part_index, attr_index) - returns double-value
  //  (attr_name, part_index, attr_index, value) - sets the new value to the attribute
  //This is slow. In the C++ code you have to get reference to
  //particles' attributes object and operate through it
  static PyObject* Bunch_partAttrValue(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 3 this is get attribute
    //if nVars == 4 this is set attribute
    int nVars = PyTuple_Size(args);

    const char* attr_name = NULL;
    int part_index = 0;
    int attr_index = 0;
    double val = 0.;

    if(nVars == 3 ||  nVars == 4){
       if(!PyArg_ParseTuple(args,"sii|d:partAttrValue",&attr_name,&part_index ,&attr_index,&val)){
          error("PyBunch - partAttrValue(attr_name,part_index,atr_index,[val]) - params. are needed");
       }
             std::string attr_name_str(attr_name);
             int bunch_size = cpp_bunch->getSize();
             int attr_size = cpp_bunch->getParticleAttributes(attr_name_str)->getAttSize();
             if(part_index >=  bunch_size || attr_size <= attr_index){
                 error("PyBunch - partAttrValue(attr_name,part_index,atr_index,[val]) - indexes out of limits! Stop!");
             }
      if(nVars == 3){
        val = cpp_bunch->getParticleAttributes(attr_name_str)->attValue(part_index,attr_index);
        return Py_BuildValue("d",val);
      }
      else{
        cpp_bunch->getParticleAttributes(attr_name_str)->attValue(part_index,attr_index) = val;
                return Py_BuildValue("d",val);
      }
    }
    else{
      error("PyBunch. You should call partAttrValue(attr_name,part_ind,attr_ind) or partAttrValue(attr_name,part_ind,attr_ind,value)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //---------------------------------------------------------------
  //
  // getSize, getSizeGlobal, getSizeGlobalFromMemory, setTotalCount
  // getCapacity
  //
  //----------------------------------------------------------------

  //returns the number of macro-particles in the bunch
  //this is implementation of the "getSize()" method
  static PyObject* Bunch_getSize(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    return Py_BuildValue("i",cpp_bunch->getSize());
  }

  //returns the number of macro-particles in the bunch in all CPUs
  //this is implementation of the "getSizeGlobal()" method
  static PyObject* Bunch_getSizeGlobal(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    return Py_BuildValue("i",cpp_bunch->getSizeGlobal());
  }

  //returns the number of macro-particles in the bunch in all CPUs
  //    that was calculated in the previous call of getSizeGlobal()
  //this is implementation of the "getSizeGlobalFromMemory()" method
  static PyObject* Bunch_getSizeGlobalFromMemory(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    return Py_BuildValue("i",cpp_bunch->getSizeGlobalFromMemory());
  }

  //returns the number of all macro-particles - alive, dead, new
  //this is implementation of the "getTotalCount()" method
  static PyObject* Bunch_getTotalCount(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    return Py_BuildValue("i",cpp_bunch->getTotalCount());
  }

  //returns the capacity of the bunch-container. It could be changed.
  //this is implementation of the "getCapacity()" method
  static PyObject* Bunch_getCapacity(PyObject *self, PyObject *Py_UNUSED(ignored)){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    return Py_BuildValue("i",cpp_bunch->getCapacity());
  }

  //---------------------------------------------------------------
  //
  // write into file or print Bunch
  //
  //----------------------------------------------------------------

  //Prints bunch into the std::cout stream
  static PyObject* Bunch_dumpBunch(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 0 dumpBunchs into std::cout
    //if nVars == 1 dumpBunchs into the file
    int nVars = PyTuple_Size(args);
    const char* file_name = NULL;
    if(nVars == 0 ||  nVars == 1){
      if(nVars == 0){
        cpp_bunch->print(std::cout);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"s:dumpBunch",&file_name)){
          error("PyBunch - dumpBunch(fileName) - a new value are needed");
        }
        cpp_bunch->print(file_name);
      }
    }
    else{
      error("PyBunch. You should call dumpBunch() or dumpBunch(file_name)");
    }

    Py_INCREF(Py_None);
    return Py_None;
  }

  //Reads bunch info from the file
  static PyObject* Bunch_readBunch(PyObject *self, PyObject *args){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    //if nVars == 1 reads all macro-particles
    //if nVars == 2 reads only specified number of macro-particles
    int nVars = PyTuple_Size(args);
    const char* file_name = NULL;
    int nParts = 0;
    if(nVars == 1 ||  nVars == 2){
      if(nVars == 1){
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"s:read",&file_name)){
          error("PyBunch - readBunch(fileName) - a file name are needed");
        }
                cpp_bunch->initBunchAttributes(file_name);
                cpp_bunch->readParticleAttributes(file_name);
        cpp_bunch->readBunchCoords(file_name);
      }
      else{
        //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
        if(!PyArg_ParseTuple(args,"si:read",&file_name,&nParts)){
          error("PyBunch - readBunch(fileName,nParts) - file name, and number of particles are needed");
        }
                cpp_bunch->initBunchAttributes(file_name);
                cpp_bunch->readParticleAttributes(file_name);
        cpp_bunch->readBunchCoords(file_name,nParts);
      }
    }
    else{
      error("PyBunch. You should call readBunch(file_name) or readBunch(file_name,nParts)");
    }
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Copy bunch attrubutes and structure to another bunch
  static PyObject* Bunch_copyEmptyBunchTo(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    PyObject* pyBunch_Target = arg;
        Bunch* cpp_target_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch_Target)->cpp_obj;
        cpp_bunch->copyEmptyBunchTo(cpp_target_bunch);
    Py_INCREF(Py_None);
    return Py_None;
  }

  //Copy bunch all info including particles coordinates and attributes to another bunch
  static PyObject* Bunch_copyBunchTo(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    PyObject* pyBunch_Target = arg;
        Bunch* cpp_target_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch_Target)->cpp_obj;
        cpp_bunch->copyBunchTo(cpp_target_bunch);
        Py_INCREF(Py_None);
    return Py_None;
  }

  //Copy particles coordinates from one bunch to another
  static PyObject* Bunch_addParticlesTo(PyObject *self, PyObject *arg){
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) self)->cpp_obj;
    PyObject* pyBunch_Target = arg;
        Bunch* cpp_target_bunch =(Bunch*) ((pyORBIT_Object *) pyBunch_Target)->cpp_obj ;
        cpp_bunch->addParticlesTo(cpp_target_bunch);
        Py_INCREF(Py_None);
    return Py_None;
  }

  //-----------------------------------------------------
  //destructor for python Bunch class
  //-----------------------------------------------------
  //this is implementation of the __del__ method
  static void Bunch_del(pyORBIT_Object* self){
        Bunch* cpp_bunch = (Bunch*) self->cpp_obj;
        delete cpp_bunch;
        self->ob_base.ob_type->tp_free((PyObject*)self);
  }

#ifdef PyORBIT_EXPERIMENTAL_WITH_NUMPY
static PyArrayObject *parse_bunch_array(PyObject *input) {
    PyArrayObject *array = (PyArrayObject *)PyArray_FROM_OTF(
            input,
            NPY_FLOAT64,
            NPY_ARRAY_IN_ARRAY
    );

    if (array == NULL) {
        return NULL;
    }

    if (PyArray_NDIM(array) != 2) {
        PyErr_SetString(PyExc_ValueError,
                        "array must be 2-dimensional with shape (nparts, 6)"
        );
        Py_DECREF(array);
        return NULL;
    }

    if (PyArray_DIM(array, 1) != 6) {
        PyErr_SetString(PyExc_ValueError,
                        "expected shape (nparts, 6) for (x, xp, y, yp, z, dE)"
        );
        Py_DECREF(array);
        return NULL;
    }

    return array; // new ref; caller MUST Py_DECREF().
}

static void append_bunch_with_PyArray(Bunch *bunch, PyArrayObject *array) {
    const npy_intp nparts = PyArray_DIM(array, 0);
    const double *data = (const double *)PyArray_DATA(array);

    for (npy_intp i = 0; i < nparts; ++i) {
        const double *coords = data + i*6;

        bunch->addParticle(coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
    }
}

PyDoc_STRVAR(
    Bunch_to_numpy_doc,
    "to_numpy($self, /, gather=False, root=0)\n"
    "--\n"
    "\n"
    "Return the particle coordinates as a NumPy array.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "gather : bool, optional\n"
    "    Collect particles from all MPI ranks. The default is False.\n"
    "root : int, optional\n"
    "    Rank that receives the collected array. The default is 0.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "numpy.ndarray or None\n"
    "    Array with shape ``(n_particles, 6)`` and dtype ``float64``.\n"
    "    When gathering, non-root ranks return None.\n"
    "\n"
    "Raises\n"
    "------\n"
    "ValueError\n"
    "    If root is not a valid MPI rank.\n"
    "OverflowError\n"
    "    If the bunch is too large for the MPI gather operation.\n"
);

static PyObject *Bunch_to_numpy(PyObject *self, PyObject *args, PyObject *kwargs) {
  static const char* kwlist[] = {"gather", "root", NULL};

  int gather = 0;
  int root = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|pi:to_numpy", kwlist, &gather, &root)) {
      return NULL;
  }

  Bunch *bunch = (Bunch *)((pyORBIT_Object *)self)->cpp_obj;

  const int rank = bunch->getMPI_Rank();
  const int size = bunch->getMPI_Size();
  const npy_intp nparts = (npy_intp)bunch->getSize();
  const npy_intp ncoords = 6;

  if (gather && (root < 0 || root >= size)) {
      PyErr_Format(PyExc_ValueError, "root must be between 0 and %d, got %d", size - 1, rank);
      return NULL;
  }

  npy_intp dims[2] = { nparts, ncoords };

  PyObject *local_array = PyArray_SimpleNew(2, dims, NPY_FLOAT64);

  // In the unlikely event that we try to gather a bunch with size larger than INT_MAX
  // or if we try to allocate a buffer to hold the bunch data and it fails on one rank
  // then propagate that result to the other ranks and raise.
  if (gather && size > 1) {
      int local_counts_ok = nparts <= INT_MAX / ncoords;
      int all_counts_ok = 0;

      ORBIT_MPI_Allreduce(&local_counts_ok, &all_counts_ok, 1, MPI_INT, MPI_MIN, bunch->getMPI_Comm_Local()->comm);

      if(!all_counts_ok) {
          PyErr_SetString(PyExc_OverflowError, "local bunch is too large for MPI_Gatherv");
          return NULL;
      }

      int local_alloc_ok = local_array != NULL;
      int all_alloc_ok = 0;

      ORBIT_MPI_Allreduce(&local_alloc_ok, &all_alloc_ok, 1, MPI_INT, MPI_MIN, bunch->getMPI_Comm_Local()->comm);

      if(!all_alloc_ok) {
          Py_XDECREF(local_array);

          if(!PyErr_Occurred()) {
             PyErr_NoMemory();
          }

          return NULL;
      }
  } else {
      if (local_array == NULL) {
          return NULL;
      }
  }

  PyArrayObject *arr_obj = (PyArrayObject*)local_array;
  double *local_data = (double*)PyArray_DATA(arr_obj);
  double **src = bunch->coordArr();

  for (npy_intp i = 0; i < nparts; ++i) {
    for (npy_intp j = 0; j < ncoords; ++j) {
      local_data[j + i*ncoords] = src[i][j];
    }
  }

  if (!gather || size == 1) {
    return local_array;
  }

#if USE_MPI > 0
  MPI_Comm comm = bunch->getMPI_Comm_Local()->comm;

  std::vector<int> global_nparts(size);
  std::vector<int> displacements(size);
  std::vector<int> recv_counts(size);

  const int local_nparts = (int)nparts;

  if (MPI_SUCCESS != MPI_Gather(&local_nparts, 1, MPI_INT, global_nparts.data(), 1, MPI_INT, root, comm)) {
      Py_DECREF(local_array);
      PyErr_SetString(PyExc_RuntimeError, "MPI_Gather failed to collect the bunch sizes across ranks");
      return NULL;
  }

  int total = 0;
  int layout_ok = 1;

  if (rank == root) {
      for (int mpi_rank = 0; mpi_rank < size; ++mpi_rank) {
          const int rank_nvalues = global_nparts[mpi_rank] * ncoords;

          if (total + rank_nvalues > INT_MAX) {
              layout_ok = 0;
              break;
          }

          displacements[mpi_rank] = total;
          recv_counts[mpi_rank] = rank_nvalues;
          total += rank_nvalues;
      }
  }

  ORBIT_MPI_Bcast(&layout_ok, 1, MPI_INT, root, comm);

  if(!layout_ok) {
      Py_DECREF(local_array);
      PyErr_SetString(PyExc_OverflowError, "global bunch is too large to collect");
      return NULL;
  }

  PyObject *global_array = NULL;

  if (rank == root) {
      npy_intp global_dims[2] = {(npy_intp)(total / ncoords), ncoords};
      global_array = PyArray_SimpleNew(2, global_dims, NPY_FLOAT64);
  }

  int global_alloc_ok = rank != root || global_array != NULL;

  ORBIT_MPI_Bcast(&global_alloc_ok, 1, MPI_INT, root, comm);

  if(!global_alloc_ok) {
      Py_DECREF(local_array);

      if(!PyErr_Occurred()) {
          PyErr_NoMemory();
      }

      return NULL;
  }

  double* global_data = rank == root ? (double *)PyArray_DATA((PyArrayObject*)global_array) : NULL;

  if(MPI_SUCCESS != MPI_Gatherv(local_data, nparts*ncoords, MPI_DOUBLE, global_data, recv_counts.data(), displacements.data(), MPI_DOUBLE, root, comm)) {
      Py_XDECREF(global_array);
      PyErr_SetString(PyExc_RuntimeError, "MPI_Gatherv failed to collect bunch.");
      return NULL;
  }

  if (rank == root) {
      return global_array;
  }

#endif // USE_MPI > 0
  Py_RETURN_NONE;
}

PyDoc_STRVAR(
    Bunch_update_from_numpy_doc,
    "update_from_numpy($self, array, /)\n"
    "--\n"
    "\n"
    "Replace the local particle coordinates from an array.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "array : array_like\n"
    "    Particle coordinates with shape ``(n_particles, 6)`` ordered as\n"
    "    ``(x, xp, y, yp, z, dE)``. Values are converted to ``float64``.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "None\n"
    "\n"
    "Raises\n"
    "------\n"
    "ValueError\n"
    "    If the array does not have shape ``(n_particles, 6)``.\n"
);

static PyObject *Bunch_update_from_numpy(PyObject *self, PyObject *arg) {
  Bunch *bunch = (Bunch *)((pyORBIT_Object *)self)->cpp_obj;

  PyArrayObject *array = parse_bunch_array(arg);

  if (array == NULL) {
      return NULL;
  }

  bunch->deleteAllParticles();
  append_bunch_with_PyArray(bunch, array);

  Py_DECREF(array);
  Py_RETURN_NONE;
}

PyDoc_STRVAR(
    Bunch_from_numpy_doc,
    "from_numpy($type, array, /)\n"
    "--\n"
    "\n"
    "Construct a bunch from particle coordinates.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "array : array_like\n"
    "    Particle coordinates with shape ``(n_particles, 6)`` ordered as\n"
    "    ``(x, xp, y, yp, z, dE)``. Values are converted to ``float64``.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "Bunch\n"
    "    A new bunch containing the supplied particles.\n"
    "\n"
    "Raises\n"
    "------\n"
    "ValueError\n"
    "    If the array does not have shape ``(n_particles, 6)``.\n"
);

static PyObject *Bunch_from_numpy(PyObject *cls, PyObject *arg) {
  PyArrayObject *array = parse_bunch_array(arg);

  if (array == NULL) {
    return NULL;
  }

  PyObject *py_bunch_obj = PyObject_CallNoArgs(cls);

  if (py_bunch_obj == NULL) {
      Py_DECREF(array);
      return NULL;
  }

  Bunch *bunch = (Bunch*)((pyORBIT_Object*)py_bunch_obj)->cpp_obj;

  append_bunch_with_PyArray(bunch, array);

  Py_DECREF(array);
  return py_bunch_obj;
}
#endif // PyORBIT_EXPERIMENTAL_WITH_NUMPY

  static PyMethodDef BunchClassMethods[] = {
    //--------------------------------------------------------
    // class Bunch wrapper                        START
    //--------------------------------------------------------
    { "getMPIComm",                     Bunch_getMPIComm                    ,METH_NOARGS,"Returns MPI Comm of this bunch"},
    { "setMPIComm",                     Bunch_setMPIComm                    ,METH_O,"Sets a new MPI Comm for this bunch"},
    { "getSyncParticle",                Bunch_getSyncParticle               ,METH_NOARGS,"Returns syncParticle class instance"},
    { "addParticle",                    Bunch_addParticle                   ,METH_VARARGS,"Adds a macro-particle to the bunch"},
    { "deleteParticle",                 Bunch_deleteParticle                ,METH_O,"Removes macro-particle from the bunch and call compress inside"},
    { "deleteParticleFast",             Bunch_deleteParticleFast            ,METH_O,"Removes macro-particle from the bunch very fast"},
    { "recoverParticle",                Bunch_recoverParticle               ,METH_O,"Recovers a particle marked for removal"},
    { "deleteAllParticles",             Bunch_deleteAllParticles            ,METH_NOARGS,"Removes all macro-particles from the bunch"},
    { "compress",                       Bunch_compress                      ,METH_NOARGS,"Compress the bunch"},
    { "x",                              Bunch_x                             ,METH_VARARGS,"Set x(index,value) or get x(index) coordinate"},
    { "y",                              Bunch_y                             ,METH_VARARGS,"Set y(index,value) or get y(index) coordinate"},
    { "z",                              Bunch_z                             ,METH_VARARGS,"Set z(index,value) or get z(index) coordinate"},
    { "px",                             Bunch_px                            ,METH_VARARGS,"Set px(index,value) or get px(index) coordinate"},
    { "py",                             Bunch_py                            ,METH_VARARGS,"Set py(index,value) or get py(index) coordinate"},
    { "pz",                             Bunch_pz                            ,METH_VARARGS,"Set pz(index,value) or get pz(index) coordinate"},
    { "dE",                             Bunch_pz                            ,METH_VARARGS,"Set dE(index,value) or get dE(index) coordinate"},
    { "xp",                             Bunch_px                            ,METH_VARARGS,"Set xp(index,value) or get xp(index) coordinate"},
    { "yp",                             Bunch_py                            ,METH_VARARGS,"Set yp(index,value) or get yp(index) coordinate"},
    { "flag",                           Bunch_flag                          ,METH_O,"Returns flag(index) for particle with index"},
    { "ringwrap",                       Bunch_ringwrap                      ,METH_O,"Perform the ring wrap. Usage: ringwrap(ring_length)"},
    { "mass",                           Bunch_mass                          ,METH_VARARGS,"Set mass(value) or get mass() the mass of particle in MeV"},
    { "classicalRadius",                Bunch_classicalRadius               ,METH_NOARGS,"Returns a classical radius of particle in [m]"},
    { "B_Rho",                          Bunch_B_Rho                         ,METH_NOARGS,"Returns B*Rho parameter of particle in [T*m]"},
    { "charge",                         Bunch_charge                        ,METH_VARARGS,"Set charge(value) or get charge() the charge of particle in e-charge"},
    { "macroSize",                      Bunch_macroSize                     ,METH_VARARGS,"Set macroSize(value) or get macroSize() the charge of particle in e-charge"},
    { "initBunchAttr",                  Bunch_initBunchAttr                 ,METH_O,"Reads and initilizes bunch attributes from a bunch file"},
    { "bunchAttrDouble",                Bunch_bunchAttrDouble               ,METH_VARARGS,"Returns and sets a double bunch attribute"},
    { "bunchAttrInt",                   Bunch_bunchAttrInt                  ,METH_VARARGS,"Returns and sets an integer bunch attribute"},
    { "bunchAttrDoubleNames",           Bunch_bunchAttrDoubleNames          ,METH_NOARGS,"Returns a list of double bunch attribute names"},
    { "bunchAttrIntNames",              Bunch_bunchAttrIntNames             ,METH_NOARGS,"Returns a list of integer bunch attribute names"},
    { "hasBunchAttrDouble",             Bunch_hasBunchAttrDouble            ,METH_O,"Returns 1 if there is a double bunch attr. with this name, 0 - otherwise"},
    { "hasBunchAttrInt",                Bunch_hasBunchAttrInt               ,METH_O,"Returns 1 if there is a int bunch attr. with this name, 0 - otherwise"},
    { "addPartAttr",                    Bunch_addPartAttr                   ,METH_VARARGS,"Adds a particles' attributes to the bunch"},
    { "removePartAttr",                 Bunch_removePartAttr                ,METH_O,"Removes a particles' attributes from the bunch"},
    { "removeAllPartAttr",              Bunch_removeAllPartAttr             ,METH_NOARGS,"Removes all particles' attributes from the bunch"},
    { "getPartAttrNames",               Bunch_getPartAttrNames              ,METH_NOARGS,"Returns all particles' attributes names in the bunch at this moment"},
    { "getPartAttrDicts",               Bunch_getPartAttrDicts              ,METH_NOARGS,"Returns dict{part. attribute name:dict{key:val}}"},
    { "getPossiblePartAttrNames",       Bunch_getPossiblePartAttrNames      ,METH_NOARGS,"Returns all possible particles' attributes names"},
    { "clearAllPartAttrAndMemorize",    Bunch_clearAllPartAttrAndMemorize   ,METH_NOARGS,"Temporary removes and memorizes all particles' attributes names"},
    { "restoreAllPartAttrFromMemory",   Bunch_restoreAllPartAttrFromMemory  ,METH_NOARGS,"Restores all particles' attributes names from memory"},
    { "hasPartAttr",                    Bunch_hasPartAttr                   ,METH_O,"Returns 1 if there is a particles' attr. with this name, 0 - otherwis"},
    { "readPartAttrNames",              Bunch_readPartAttrNames             ,METH_O,"Returns a tuple with particles' attr. names in the bunch file"},
    { "readPartAttrDicts",              Bunch_readPartAttrDicts             ,METH_O,"Returns a dict{attr_name:dicts{param_name:val}} in the bunch file"},
    { "readPartAttr",                   Bunch_readPartAttr                  ,METH_O,"Initializes the particles' attr. from the bunch file"},
    { "getPartAttrSize",                Bunch_getPartAttrSize               ,METH_O,"Returns the number of variables in the particles' attributes with a particular name"},
    { "partAttrValue",                  Bunch_partAttrValue                 ,METH_VARARGS,"Sets or returns a particles' attribute value"},
    { "getSize",                        Bunch_getSize                       ,METH_NOARGS,"Returns number of macro-particles"},
    { "getSizeGlobal",                  Bunch_getSizeGlobal                 ,METH_NOARGS,"Returns number of macro-particles in all CPUs"},
    { "getSizeGlobalFromMemory",        Bunch_getSizeGlobalFromMemory       ,METH_NOARGS,"Returns number of macro-particles in all CPUs from memory"},
    { "getTotalCount",                  Bunch_getTotalCount                 ,METH_NOARGS,"Returns number of all particles - alive,dead,new"},
    { "getCapacity",                    Bunch_getCapacity                   ,METH_NOARGS,"Returns the capacity of the bunch-contaiter"},
    { "dumpBunch",                      Bunch_dumpBunch                     ,METH_VARARGS,"Prints the bunch info into a standart output stream or file"},
    { "readBunch",                      Bunch_readBunch                     ,METH_VARARGS,"Reads the bunch info from a file"},
    { "copyEmptyBunchTo",               Bunch_copyEmptyBunchTo              ,METH_O,"Copy bunch attrubutes and structure to another bunch"},
    { "copyBunchTo",                    Bunch_copyBunchTo                   ,METH_O,"Copy bunch all info including particles coordinates and attributes to another bunch"},
    { "addParticlesTo",                 Bunch_addParticlesTo                ,METH_O,"Copy particles coordinates from one bunch to another"},
#ifdef PyORBIT_EXPERIMENTAL_WITH_NUMPY
    { "to_numpy",                       _PyCFunction_CAST(Bunch_to_numpy)                      ,METH_VARARGS | METH_KEYWORDS, Bunch_to_numpy_doc },
    { "update_from_numpy",              Bunch_update_from_numpy             ,METH_O, Bunch_update_from_numpy_doc },
    { "from_numpy",                     Bunch_from_numpy                    ,METH_O | METH_CLASS, Bunch_from_numpy_doc },
#endif
    {NULL,NULL}
    //--------------------------------------------------------
    // class Bunch wrapper                        STOP
    //--------------------------------------------------------
  };

    // defenition of the memebers of the python Bunch wrapper class
    // they will be vailable from python level
    static PyMemberDef BunchClassMembers [] = {
        {NULL}
    };

    //new python Bunch wrapper type definition
    static PyTypeObject pyORBIT_Bunch_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "Bunch", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) Bunch_del , /*tp_dealloc*/
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
        "The Bunch python wrapper", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        BunchClassMethods, /* tp_methods */
        BunchClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) Bunch_init, /* tp_init */
        0, /* tp_alloc */
        Bunch_new, /* tp_new */
    };

  static PyMethodDef BunchModuleMethods[] = { {NULL,NULL} };


#ifdef __cplusplus
extern "C" {
#endif

  static struct PyModuleDef cModPyDem =
  {
      PyModuleDef_HEAD_INIT,
      "bunch", "Bunch class",
      -1,
      BunchModuleMethods
  };

  /* The name of the function was changed to avoid collision with PyImport magic naming */
  PyMODINIT_FUNC initbunch(void) {
  #ifdef PyORBIT_EXPERIMENTAL_WITH_NUMPY
    if (ensure_numpy() != 0) {
      return NULL;
    }
  #endif // PyORBIT_EXPERIMENTAL_WITH_NUMPY
      //check that the Bunch wrapper is ready
      if(PyType_Ready(&pyORBIT_Bunch_Type) < 0) return NULL;
      Py_INCREF(&pyORBIT_Bunch_Type);
      PyObject* module = PyModule_Create(&cModPyDem);
      PyModule_AddObject(module, "Bunch", (PyObject *)&pyORBIT_Bunch_Type);
      //add the SyncParticle python class
      wrap_orbit_syncpart::initsyncpart(module);
      wrap_bunch_twiss_analysis::initbunchtwissanalysis(module);
      wrap_bunch_tune_analysis::initbunchtuneanalysis(module);
      wrap_synch_part_redefinition::initsynchpartredefinition(module);
      return module;
  }

	PyObject* getBunchType(const char* name){
		PyObject* mod = PyImport_ImportModule("orbit.core.bunch");
		PyObject* pyType = PyObject_GetAttrString(mod,name);
		Py_DECREF(mod);
		Py_DECREF(pyType);
		return pyType;
	}


#ifdef __cplusplus
}
#endif

//end of namespace wrap_orbit_bunch
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
