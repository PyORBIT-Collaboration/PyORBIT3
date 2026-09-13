#include "main/pyORBIT_Object.hh"

#include <unordered_map>

namespace {

std::unordered_map<void*, PyObject*> wrappers;

}

namespace pyorbit {

void registerPyWrapper(void* cpp_obj, PyObject* py_wrapper)
{
  wrappers[cpp_obj] = py_wrapper;
}

void unregisterPyWrapper(void* cpp_obj, PyObject* py_wrapper)
{
  auto wrapper = wrappers.find(cpp_obj);
  if(wrapper != wrappers.end() && wrapper->second == py_wrapper) {
    wrappers.erase(wrapper);
  }
}

PyObject* getPyWrapper(void* cpp_obj)
{
  auto wrapper = wrappers.find(cpp_obj);
  return wrapper == wrappers.end() ? nullptr : wrapper->second;
}

}
