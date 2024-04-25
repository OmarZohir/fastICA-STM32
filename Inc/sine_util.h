#ifndef __SINE_UTIL_H__
#define __SINE_UTIL_H__

#include "arm_math.h"
#include "stdlib.h"
#include "stm32l475e_iot01_qspi.h"
#include <cmath>

#define QSPI_FLASH_BLOCK_SIZE (64*1024)
#define MIXER_BLOCK 100

/* Generate floating point samples of a sine wave at a given sampling rate and frequency. Stores in memory. */
int sine_to_buffer(int frequency, int sampling_rate, int duration_ms, int block_number);

int sine_to_buffer_sram(int frequency, int sampling_rate, int duration_ms, float * buffer);

/* Convert a normalized (-1 to 1) floating point and return a 12-bit right aligned value for the DAC. */
uint32_t convert_f32_to_u12(float value);

void mix_from_flash(float * srcAddr1, float * srcAddr2, uint32_t size, uint32_t destAddressFlash);

/* Generate an already-mixed already normalized sine of two samples. */
uint32_t generate_mixed_sines(float f1, float f2, int sampling_rate, int duration_ms, uint32_t pDest);

#endif
