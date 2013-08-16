#include "jip_dispatcher.h"
#include "Python.h"


bool _update_files(PyObject* list, FILE** target, uint64_t elements){
  PyObject* item;
  FILE* current;
  uint64_t i = 0;
  for (i = 0; i < elements; i++) {
    item = PySequence_Fast_GET_ITEM(list, i);    
    if(PyString_Check(item)){
      current = fopen(PyString_AsString(item), i==0 ? "r": "w");
      if(current == NULL){
        fprintf(stderr, "ERROR WHILE OPENING FILE!\n");
        // unable to open file
        PyErr_SetFromErrno(PyExc_OSError);
        return false;
      }
    }else{
      current = PyFile_AsFile(item);
    }
    target[i] = current;
  }
  return true;
}

static PyObject* dispatch_streams(PyObject *self, PyObject *args){
  /* Parse argument tuple and open FILE* for input and
   * output. This should raise an exception if something fails
   */  
  int i = 0;
  PyObject* seq = PySequence_Fast(args, "Expexted the argument list");
  uint64_t num_sources = 0;
  uint64_t num_targets = 0;
  uint64_t len = PySequence_Size(args);
  if(seq == NULL || len != 2){
    PyErr_SetString(PyExc_ValueError, "Expected exactly two arguments!");
    Py_DECREF(seq);
    return NULL;
  }
  PyObject* py_sources = PySequence_Fast_GET_ITEM(seq, 0);
  if(!PySequence_Check(py_sources)){
    PyErr_SetString(PyExc_ValueError, "Argument 0 is not a list!");
    Py_DECREF(seq);
    return NULL;
  }
  PyObject* py_targets = PySequence_Fast_GET_ITEM(seq, 1);
  if(!PySequence_Check(py_targets)){
    PyErr_SetString(PyExc_ValueError, "Argument 1 is not a list!");
    Py_DECREF(seq);
    return NULL;
  }
  num_sources = PySequence_Size(py_sources);
  num_targets = PySequence_Size(py_targets);
  FILE** source_f = malloc(num_sources * sizeof(FILE*));
  FILE** target_f = malloc(num_targets * sizeof(FILE*));
  
  bool error = false;
  if(!_update_files(py_sources, source_f, num_sources)){
    error = true;
  }
  if(!_update_files(py_targets, target_f, num_targets)){
    error = true;
  }

  if(!error){
    dispatch(source_f, target_f, num_sources, num_targets);
  }


  //for (i=0;i<num_sources;i++) {
  //  fclose(source_f[i]);
  //}
  //for (i=0;i<num_targets;i++) {
  //  fclose(target_f[i]);
 // }


  // cleanup and free the file lists
  free(target_f);
  free(source_f);
  Py_DECREF(seq);

  if(error) return NULL;
  Py_RETURN_NONE;
}

static PyMethodDef DispatchMethods[] = {
    {"dispatch",  dispatch_streams, METH_VARARGS, 
     "Dispatch data from source to all targets. The method "
     "excepts strings and open file handles you hev to specify at least "
     "one source and one target."},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initdispatcher(void){
      (void) Py_InitModule("dispatcher", DispatchMethods);
}
