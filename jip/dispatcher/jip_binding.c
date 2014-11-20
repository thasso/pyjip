#include "Python.h"
#include "jip_dispatcher.h"

#if PY_MAJOR_VERSION >= 3
#define PyString_Check PyUnicode_Check
#define PyString_AsString PyUnicode_AS_DATA
#endif


bool _update_files(PyObject* list, FILE** target, uint64_t elements){
  PyObject* item;
  FILE* current;
  uint64_t i = 0;
  int fd1 ;
  PyObject *mode_obj1 ;
  PyObject *mode_byte_obj1;
  char *mode1;
  for (i = 0; i < elements; i++) {
    item = PySequence_Fast_GET_ITEM(list, i);
    if(Py_None == item){
      target[i] = NULL;
    }else{
      if(PyString_Check(item)){
        current = fopen(PyString_AsString(item), i==0 ? "r": "w");
        if(current == NULL){
          fprintf(stderr, "ERROR WHILE OPENING FILE!\n");
          // unable to open file
          PyErr_SetFromErrno(PyExc_OSError);
          return false;
        }
      }else{
#if PY_MAJOR_VERSION >= 3
        fd1 = PyObject_AsFileDescriptor(item);
        mode_obj1 = PyObject_GetAttrString(item, "mode");
        mode_byte_obj1 = PyUnicode_AsUTF8String(mode_obj1);
        mode1 = PyBytes_AsString(mode_byte_obj1);
        current = fdopen(fd1, mode1);
        Py_XDECREF(mode_obj1);
        Py_XDECREF(mode_byte_obj1);
#else
        current = PyFile_AsFile(item);
#endif
      }
      target[i] = current;
    }
  }
  return true;
}

PyObject* dispatch_streams_(PyObject *self, PyObject *args, int (*fun)(FILE**, FILE **, FILE **, uint64_t)){
  /* Parse argument tuple and open FILE* for input and
   * output. This should raise an exception if something fails
   */
  PyObject* seq = PySequence_Fast(args, "Expected the argument list");
  PyObject* py_sources = PySequence_Fast_GET_ITEM(seq, 0);
  PyObject* py_targets_1 = PySequence_Fast_GET_ITEM(seq, 1);
  PyObject* py_targets_2 = PySequence_Fast_GET_ITEM(seq, 2);
  long num_sources = 0;
  long num_targets_1 = 0;
  long num_targets_2 = 0;
  uint64_t len = PySequence_Size(args);
  FILE** source_f;
  FILE** target_f_1;
  FILE** target_f_2;
  bool error;

  if(seq == NULL || len != 3){
    PyErr_SetString(PyExc_ValueError, "Expected exactly three arguments!");
    Py_DECREF(seq);
    return NULL;
  }
  if(!PySequence_Check(py_sources)){
    PyErr_SetString(PyExc_ValueError, "Argument 0 is not a list!");
    Py_DECREF(seq);
    return NULL;
  }
  if(!PySequence_Check(py_targets_1)){
    PyErr_SetString(PyExc_ValueError, "Argument 1 is not a list!");
    Py_DECREF(seq);
    return NULL;
  }
  if(!PySequence_Check(py_targets_2)){
    PyErr_SetString(PyExc_ValueError, "Argument 2 is not a list!");
    Py_DECREF(seq);
    return NULL;
  }
  num_sources = PySequence_Size(py_sources);
  num_targets_1 = PySequence_Size(py_targets_1);
  num_targets_2 = PySequence_Size(py_targets_2);

  source_f = malloc(num_sources * sizeof(FILE*));
  target_f_1 = malloc(num_targets_1 * sizeof(FILE*));
  target_f_2 = malloc(num_targets_2 * sizeof(FILE*));

  error = false;
  if(!_update_files(py_sources, source_f, num_sources)){
    error = true;
  }
  if(!_update_files(py_targets_1, target_f_1, num_targets_1)){
    error = true;
  }
  if(!_update_files(py_targets_2, target_f_2, num_targets_2)){
    error = true;
  }

  if(!error){
    // for a child process for the dispatcher
    pid_t childPID = fork();
    if(childPID >= 0){
      if(childPID==0){
        // child
        fun(source_f, target_f_1, target_f_2, num_sources);
        free(target_f_1);
        free(target_f_2);
        free(source_f);
        exit(0);
      }
    }else{
      free(target_f_1);
      free(target_f_2);
      free(source_f);
      PyErr_SetString(PyExc_ValueError, "Unable to create dispatcher process");
      error = true;
    }

  }
  // cleanup
  Py_DECREF(seq);
  if(error) return NULL;
  Py_RETURN_NONE;
}

static PyObject* dispatch_streams(PyObject *self, PyObject *args){
  return dispatch_streams_(self, args, &dispatch);
}

static PyObject* dispatch_streams_fanout(PyObject *self, PyObject *args){
  return dispatch_streams_(self, args, &dispatch_fanout);
}

static PyObject* dispatch_streams_fanin(PyObject *self, PyObject *args){
  return dispatch_streams_(self, args, &dispatch_fanin);
}

static PyMethodDef DispatchMethods[] = {
    {"dispatch",  dispatch_streams, METH_VARARGS,
     "Dispatch data from source to all targets. The method "
     "excepts strings and open file handles you have to specify at least "
     "one source and one target."},
    {"dispatch_fanout",  dispatch_streams_fanout, METH_VARARGS,
     "Dispatch data from source to all targets. The method "
     "excepts strings and open file handles you have to specify at least "
     "one source and one target."},
    {"dispatch_fanin",  dispatch_streams_fanin, METH_VARARGS,
     "Dispatch data from source to all targets. The method "
     "excepts strings and open file handles you have to specify at least "
     "one source and one target."},
    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "dispatcher",        /* m_name */
    "JIP Dispatcher",    /* m_doc */
    -1,                  /* m_size */
    DispatchMethods,     /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_dispatcher(void) {
    return PyModule_Create(&moduledef);
}
#else
PyMODINIT_FUNC initdispatcher(void){
      (void) Py_InitModule("dispatcher", DispatchMethods);
}
#endif
