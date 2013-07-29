#include "jip_dispatcher.h"
#include "Python.h"

static PyObject* dispatch_streams(PyObject *self, PyObject *args){
  /* Parse argument tuple and open FILE* for input and
   * output. This should raise an exception if something failes
   */  
  PyObject* item;
  FILE* current;
  int i, assigned = 0;
  PyObject* seq = PySequence_Fast(args, "Expexted the argument list");
  uint64_t len = PySequence_Size(args);
  if(len < 2){
    PyErr_SetString(PyExc_ValueError, "Expected at least two arguments!");
    Py_DECREF(seq);
    return NULL;
  }
  FILE* source_f = NULL;
  FILE** target_f = malloc((len-1) * sizeof(FILE*));
  bool error = false;
  for (i = 0; i < len; i++) {
    item = PySequence_Fast_GET_ITEM(seq, i);    
    if(PyString_Check(item)){
      current = fopen(PyString_AsString(item), i==0 ? "r": "w");
      if(current == NULL){
        // unable to open file
        PyErr_SetFromErrno(PyExc_OSError);
        error = true;
        break;
      }
    }else{
      current = PyFile_AsFile(item);
    }
    if(i==0){
      source_f = current;
    }else{
      target_f[i-1] = current;
      assigned++;
    }
    /* DON'T DECREF item here */
  }
  if(!error){
    dispatch(source_f, target_f, (len-1));
  }else{
    // close streams
    if(source_f != NULL){
      fclose(source_f);
    }
    for (i = 0; i < assigned; i++) {
      fclose(target_f[i]);
    }
  }
  // free targets
  free(target_f);
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
