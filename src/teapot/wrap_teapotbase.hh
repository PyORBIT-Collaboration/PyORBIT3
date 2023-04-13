#ifdef __cplusplus
extern "C"
{
#endif

namespace wrap_teapotbase
{
    PyMODINIT_FUNC initteapotbase(void);
    PyObject* getBaseTEAPOTType(char* name);
}

#ifdef __cplusplus
}
#endif

