#ifdef __cplusplus
extern "C"
{
#endif

namespace wrap_errorbase
{
  PyMODINIT_FUNC initerrorbase(void);
  PyObject* getBaseERRORType(char* name);
}

#ifdef __cplusplus
}
#endif
