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
int dispatch(FILE** source, FILE** targets_1, FILE**  targets_2, const uint64_t num_elements);
int dispatch_fanout(FILE** source, FILE** targets_1, FILE**  targets_2, const uint64_t num_elements);
int dispatch_fanin(FILE** source, FILE** targets_1, FILE**  targets_2, const uint64_t num_elements);

#endif /* JIP_DISPATCHER_H_ */
