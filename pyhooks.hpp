#pragma once

#include "c_common.hpp"

struct _object;
typedef _object PyObject;
typedef DiffeqFuncType DiffeqMethod;
typedef EventFunc PyEventMethod;

inline bool import_CyRK__cy__pysolver_cyhook()
{
    return true;
}

inline int call_diffeq_from_cython(PyObject* x, DiffeqMethod y)
{
    return 1;
}

inline void Py_XINCREF(PyObject* x)
{
}

inline void Py_XDECREF(PyObject* x)
{
}

inline int call_pyevent_from_cython(PyObject* py_instance, PyEventMethod pyevent_method, size_t event_index, double t, double* y_ptr)
{
    return 1;
}
