#include "sine_util.h"

/*float * sine_to_buffer(int frequency, int sampling_rate, int duration_ms)
{
	int sample_number;
	size_t required_buffer = sampling_rate * duration_ms / 1000 * sizeof(float);
	float * buf = (float *) malloc(required_buffer);
	if (buf == NULL) {
		return NULL;
	}
	float sin_constant = 2000 * PI * frequency / (duration_ms);
	for (sample_number = 0; sample_number < sampling_rate * duration_ms / 1000; sample_number++) {
		buf[sample_number] = arm_sin_f32(sin_constant * sample_number);
	}
	return buf;
}*/
int sine_to_buffer(int frequency, int sampling_rate, int duration_ms, int block_number)
{
	size_t sample_number;
	float omega = 2 * PI * frequency;
	float sr = (float) sampling_rate;
	size_t total_samples = sampling_rate * duration_ms / 1000;
	for (sample_number = 0; sample_number < total_samples; sample_number++) {
		float temp = arm_sin_f32(omega * sample_number / sr);
		if (BSP_QSPI_Write((uint8_t*)&temp, (QSPI_FLASH_BLOCK_SIZE * block_number) +  sample_number * sizeof(temp), sizeof(temp)) != QSPI_OK) {
			_Error_Handler(__FILE__, __LINE__);
		}
	}
	return sample_number;
}

int sine_to_buffer_sram(int frequency, int sampling_rate, int duration_ms, float * buffer) {
	int sample_number;
	float omega = 2 * PI * frequency;
	float sr = (float) sampling_rate;
	int total_samples = sampling_rate * duration_ms / 1000;
	for (sample_number = 0; sample_number < total_samples; sample_number++) {
		buffer[sample_number] = arm_sin_f32(omega * (float)sample_number / sr);
	}
	return sample_number;
}

uint32_t convert_f32_to_u12(float value)
{
	return (uint32_t) ((value + 1)/2.0f * 2095.0f); // 2^12 / 2
}

uint32_t generate_mixed_sines(float f1, float f2, int sampling_rate, int duration_ms, uint32_t pDest)
{
	const float mixer[4] = {4, 3, 1, 2};
	
	size_t sample_number;
	float omega1 = 2 * PI * f1;
	float omega2 = 2 * PI * f2;
	size_t total_samples = sampling_rate * duration_ms / 1000;
	float sr = (float) sampling_rate;
	float temp1, temp2, res1, res2;
	for (sample_number = 0; sample_number < total_samples; sample_number++) {
		temp1 = arm_sin_f32(omega1 * sample_number / sr);
		temp2 = arm_sin_f32(omega2 * sample_number / sr);
		res1 = mixer[0] * temp1 + mixer[1] * temp2;
		res2 = mixer[2] * temp1 + mixer[3] * temp2;
		res1 = res1 / (mixer[0] + mixer[1]);
		res2 = res2 / (mixer[2] + mixer[3]);
		if (BSP_QSPI_Write((uint8_t*)&res1, pDest + (sample_number * sizeof(res1)), sizeof(res1)) != QSPI_OK) {
			_Error_Handler(__FILE__, __LINE__);
		}
		if (BSP_QSPI_Write((uint8_t*)&res2, pDest + total_samples * sizeof(res2) + (sample_number * sizeof(res2)), sizeof(res2)) != QSPI_OK) {
			_Error_Handler(__FILE__, __LINE__);
		}
	}
	return total_samples;
}


void mix_from_flash(float * srcAddr1, float * srcAddr2, uint32_t entries, uint32_t destAddressFlash)
{
	float mixer[4] = {4, 3, 1, 2};	
	uint32_t counter = 0;
	for (counter = 0; counter < entries; counter++) {
		float s1, s2;
		BSP_QSPI_Read((uint8_t *)&s1, (uint32_t)srcAddr1 + (counter * sizeof(float)), sizeof(float));
		BSP_QSPI_Read((uint8_t *)&s2, (uint32_t)srcAddr2 + (counter * sizeof(float)), sizeof(float));
		float res1 = mixer[0] * s1 + mixer[1] * s2;
		float res2 = mixer[2] * s1 + mixer[3] * s2;
		res1 = res1 / (mixer[0] + mixer[1]);
		res2 = res2 / (mixer[2] + mixer[3]);
		
		BSP_QSPI_Write((uint8_t *)&res1, destAddressFlash + (counter * sizeof(float)), sizeof(float));
		BSP_QSPI_Write((uint8_t *)&res2, destAddressFlash + ((counter + entries) * sizeof(float)), sizeof(float));
	}	
}
