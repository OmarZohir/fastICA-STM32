#ifndef __FAST_ICA_H__
#define __FAST_ICA_H__
#include "stdlib.h"
#include "arm_math.h"
#include "stm32l475e_iot01_qspi.h"

#define MAX_ITERATIONS 100
#define STATUS

extern const size_t processing_chunk_size;
extern const float EPSILON;

/* Prewhiten data in FLASH. Takes flash address and the size of the data in BYTES.
   Data is given in format where there are two rows where the first row appears entirely before
	 the second. */
void prewhiten_f32(uint32_t flash_address, uint32_t pDest, size_t size,float* whitening_data);

/* Perform FastICA on 2xn matrix in FLASH. */
void fast_ica_f32(uint32_t pSrc, uint32_t pDest, size_t size);

#endif
