#include "jip_dispatcher.h"

void dispatch(FILE** source, FILE** target, 
              const uint64_t num_sources,
              const uint64_t num_targets){
  char* buffer = malloc(JIP_BUFFER_SIZE * sizeof(char));
  size_t bytes = 0;
  uint64_t i = 0;
  uint64_t j = 0;
  for(j=0; j<num_sources;j++){
    while((bytes = fread(buffer, 1, JIP_BUFFER_SIZE, source[j])) != 0){
      for(i=0;i<num_targets; i++){
        if(fwrite(buffer, 1, bytes, target[i]) != bytes){
          // write error
        }
      }  
    }
  }
  free(buffer);
}
