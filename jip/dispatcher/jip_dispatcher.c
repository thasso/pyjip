#include "jip_dispatcher.h"

int dispatch(FILE** sources, FILE** targets_1, FILE** targets_2,
              const uint64_t num_elements){
  char* buffer = malloc(JIP_BUFFER_SIZE * sizeof(char));
  size_t bytes = 0;
  uint64_t i = 0;
  uint64_t num_done = 0;
  int result = 0;
  while(true){
    if(num_done==num_elements) break;
    i++;
    if(i >= num_elements) i=0;
    // read from source i
    if(sources[i] != NULL){
      bytes = fread(buffer, 1, JIP_BUFFER_SIZE, sources[i]);
      if(bytes == 0){
        num_done++;
        sources[i] = NULL;
      }else{
        if(targets_1[i] != NULL){
          if(fwrite(buffer, 1, bytes, targets_1[i]) != bytes){
            fprintf(stderr, "Error while writing to target!\n");
            result=1;
            break;
          }
        }
        if(targets_2[i] != NULL){
          if(fwrite(buffer, 1, bytes, targets_2[i]) != bytes){
            fprintf(stderr, "Error while writing to target!\n");
            result=1;
            break;
          }
        }
      }
    }
  }
  free(buffer);
  return result;
}


int dispatch_fanout(FILE** sources, FILE** targets_1, FILE** targets_2,
                     const uint64_t num_elements){
  char* buffer = malloc(JIP_BUFFER_SIZE * sizeof(char));
  size_t bytes = 0;
  uint64_t j = 0;
  int result = 0;
  size_t written = 0;
  while((bytes = fread(buffer, 1, JIP_BUFFER_SIZE, sources[0])) != 0){
    written = 0;
    for(j=0; j<num_elements; j++){      
      if(targets_1[j] != NULL){
        if(fwrite(buffer, 1, bytes, targets_1[j]) != bytes){
          fprintf(stderr, "Error while writing to target!\n");
          result=1;
          break;
        }else{
          written += bytes;
        }
      }
      if(targets_2[j] != NULL){
        if(fwrite(buffer, 1, bytes, targets_2[j]) != bytes){
          fprintf(stderr, "Error while writing to target!\n");
          result=1;
          break;
        }else{
          written += bytes;
        }
      }      
    }
    if(written == 0){
      result=1;
    }
    if(result==1)break;
  }
  free(buffer);
  return result;
}

int dispatch_fanin(FILE** sources, FILE** targets_1, FILE** targets_2,
                   const uint64_t num_elements){
  char* buffer = malloc(JIP_BUFFER_SIZE * sizeof(char));
  size_t bytes = 0;
  uint64_t i = 0;
  int result = 0;

  for(i=0; i<num_elements; i++){
    if(sources[i]==NULL) continue;
    while((bytes = fread(buffer, 1, JIP_BUFFER_SIZE, sources[i])) != 0){
      // write to target_1 0 -> the fan in target
      if(targets_1[0] != NULL){
        if(fwrite(buffer, 1, bytes, targets_1[0]) != bytes){
          fprintf(stderr, "Error while writing to target!\n");
          result=1;
          break;
        }
      }
      // write to the i'th direct output
      if(targets_2[i] != NULL){
        if(fwrite(buffer, 1, bytes, targets_2[i]) != bytes){
          fprintf(stderr, "Error while writing to target!\n");
          result=1;
          break;
        }
      }
    }
    if(result == 1) break;
  }
  free(buffer);
  return result;
}


