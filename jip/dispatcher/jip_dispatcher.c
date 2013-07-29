#include "jip_dispatcher.h"

void dispatch(FILE* source, FILE** target, 
              const uint64_t num_targets){
  char* buffer = malloc(JIP_BUFFER_SIZE * sizeof(char));
  size_t bytes = 0;
  uint64_t i = 0;

  while((bytes = fread(buffer, 1, JIP_BUFFER_SIZE, source)) != 0){
    for(i=0;i<num_targets; i++){
      if(fwrite(buffer, 1, bytes, target[i]) != bytes){
        // write error
      }
    }  
  }
  // close source
  fclose(source);
  // close all targets
  for(i=0;i<num_targets;i++) fclose(target[i]);
  free(buffer);
}

int main(){
  FILE** target = malloc(2 * sizeof(FILE*));
  target[0] = stdout;
  target[1] = stderr;
  dispatch(stdin, target,2);
}

