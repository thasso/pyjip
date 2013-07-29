#ifndef JIP_DISPATCHER_H_
#define JIP_DISPATCHER_H_

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <float.h>
#include <inttypes.h>
#include <stdlib.h>

#define JIP_BUFFER_SIZE 4096

/**
 * Read bytes from source file and dispatch the read bytes to
 * al the targets
 */
void dispatch(FILE* source, FILE** target, const uint64_t num_targets);

#endif /* JIP_DISPATCHER_H_ */
