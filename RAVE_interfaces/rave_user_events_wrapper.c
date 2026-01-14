#include <Python.h>
#include "rave_user_events.h"

static PyObject* py_rave_event_and_value(PyObject* self, PyObject* args) {
    int e,v;
    if (!PyArg_ParseTuple(args, "ii", &e, &v)) return NULL;
    rave_event_and_value(e,v); // Call the C function
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_restart_trace(PyObject* self, PyObject* args) {
    rave_restart_trace(); // Call the C function
    Py_RETURN_NONE;  // Return None in Python
}
//Maintaining old start/stop functions instead of enable/disable
static PyObject* py_rave_start_trace(PyObject* self, PyObject* args) {
    rave_enable_trace(); // Call the C function
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_stop_trace(PyObject* self, PyObject* args) {
    rave_disable_trace(); // Call the C function
    Py_RETURN_NONE;  // Return None in Python
}
//New naming convention:
static PyObject* py_rave_enable_trace(PyObject* self, PyObject* args) {
    rave_enable_trace(); // Call the C function
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_disable_trace(PyObject* self, PyObject* args) {
    rave_disable_trace(); // Call the C function
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_enable_regions(PyObject* self, PyObject* args) {
    rave_enable_regions(); // Call the C function
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_disable_regions(PyObject* self, PyObject* args) {
    rave_disable_regions(); // Call the C function
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_enable(PyObject* self, PyObject* args) {
    rave_enable()
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_disable(PyObject* self, PyObject* args) {
    rave_disable()
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_name_event(PyObject* self, PyObject* args) {
    int e;
		char * n;
    if (!PyArg_ParseTuple(args, "is", &e, &n)) return NULL;
		rave_name_event(e,n);
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_name_value(PyObject* self, PyObject* args) {
    int e,v;
		char * n;
    if (!PyArg_ParseTuple(args, "iis", &e, &v, &n)) return NULL;
		rave_name_value(e,v,n);
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_begin_region(PyObject* self, PyObject* args) {
		const char * name;
    if (!PyArg_ParseTuple(args, "s", &name)) return NULL;
		rave_begin_region(name); //Call the C function
    Py_RETURN_NONE;  // Return None in Python
}
static PyObject* py_rave_end_region(PyObject* self, PyObject* args) {
		const char * name;
    if (!PyArg_ParseTuple(args, "s", &name)) return NULL;
		rave_end_region(name); //Call the C function
    Py_RETURN_NONE;  // Return None in Python
}


// Define methods in the module
static PyMethodDef RaveUserEventsMethods[] = {
    {"rave_event_and_value", py_rave_event_and_value, METH_VARARGS, ""},
    {"rave_restart_trace", py_rave_restart_trace, METH_VARARGS, ""},
    {"rave_start_trace", py_rave_start_trace, METH_VARARGS, ""},
    {"rave_stop_trace", py_rave_stop_trace, METH_VARARGS, ""},
    {"rave_begin_region", py_rave_stop_trace, METH_VARARGS, ""},
    {"rave_end_region", py_rave_stop_trace, METH_VARARGS, ""},
    {"rave_name_event", py_rave_name_event, METH_VARARGS, ""},
    {"rave_name_value", py_rave_name_value, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}  // Sentinel to indicate the end of the array
};

// Module definition structure
static struct PyModuleDef rave_user_events_module = {
    PyModuleDef_HEAD_INIT,
    "rave_user_events",  // Module name
    "Python interface for the sdv_trace C library",
    -1,
    RaveUserEventsMethods
};

// Module initialization function
PyMODINIT_FUNC PyInit_rave_user_events(void) {
    return PyModule_Create(&rave_user_events_module);
}

