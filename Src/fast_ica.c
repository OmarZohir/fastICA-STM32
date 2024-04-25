#include "fast_ica.h"

const size_t processing_chunk_size = 100 * sizeof(float);
const float EPSILON = 0.0001;

void prewhiten_f32(uint32_t flash_location,  uint32_t pDest, size_t size, float* whitening_data)
{
	#ifdef STATUS
	printf("\tPrewhitening...");
	#endif
	// Store size of row
	const size_t row_size = size / 2;
	// Get FLASH address of second row
	uint32_t second_location = flash_location + row_size;
	
	/* CENTER THE ROWS */
	// Get sums of rows.
	float first_sum = 0.0f, second_sum = 0.0f, temp;
	size_t index;
	for (index = 0; index < row_size; index += sizeof(float)) {
		BSP_QSPI_Read((uint8_t*)&temp, second_location + index, sizeof(float));
		second_sum += temp;
		BSP_QSPI_Read((uint8_t*)&temp, flash_location + index, sizeof(float));
		first_sum += temp;
	}
	// Normalize sums
	float entries = (float)(row_size / sizeof(float));
	first_sum /= entries;
	second_sum /= entries;
	

	// Center rows, begin calculating covariance
	float cov_data[4] = { 0, 0, 0, 0 };	// Only 3 entries calculated directly because entry (1,2) and (2,1) are equivalent.
	float temp2;
	for (index = 0; index < row_size; index += sizeof(float)) {
		BSP_QSPI_Read((uint8_t*)&temp, flash_location + index, sizeof(float));
		temp -= first_sum;
		cov_data[0] += temp * temp;
		BSP_QSPI_Write((uint8_t*)&temp, pDest + index, sizeof(float));
		BSP_QSPI_Read((uint8_t*)&temp2, second_location + index, sizeof(float));
		temp2 -= second_sum;
		cov_data[1] += temp * temp2;
		cov_data[3] += temp2 * temp2;
		BSP_QSPI_Write((uint8_t*)&temp2, pDest + row_size + index, sizeof(float));
	}
	// NOTE done with source data
	// Finish calculating covariance matrix
	cov_data[0] /= (entries - 1);
	cov_data[1] /= (entries - 1);
	cov_data[3] /= (entries - 1);
	cov_data[2] = cov_data[1];

#ifdef DEBUG	
	printf("covariance matrix = \n%f %f\n%f %f\n\n",cov_data[0],cov_data[1],cov_data[2],cov_data[3]);
#endif
	uint32_t pDest2=pDest+size;
	
	// Calculate eigenvalues and vectors
	float tr = cov_data[0] + cov_data[3];
	float det = cov_data[0] * cov_data[3] - cov_data[1] * cov_data[2];
	arm_sqrt_f32(tr * tr - 4 * det, &temp);
	float eigenval_data[4] = {(tr - temp) / 2.0f, 0, 0, (tr + temp) / 2.0f };
	float eigenvec_data[4] = {//[C D-a][D-a;-C] = [0]
		cov_data[3] - eigenval_data[0],
		cov_data[3] - eigenval_data[3],
		-cov_data[2],
		-cov_data[2]
	};
	// Normalize eigenvectors
	arm_sqrt_f32(eigenvec_data[0]*eigenvec_data[0] + eigenvec_data[2] * eigenvec_data[2], &temp);
	eigenvec_data[0] /= temp;
	eigenvec_data[2] /= temp;
	arm_sqrt_f32(eigenvec_data[1] * eigenvec_data[1] + eigenvec_data[3] * eigenvec_data[3], &temp);
	eigenvec_data[1] /= temp;
	eigenvec_data[3] /= temp;
	
	// Get sqrt of eigenvalue matrix since we only use the sqrt from here on.
	arm_sqrt_f32(eigenval_data[0], &(eigenval_data[0]));
	arm_sqrt_f32(eigenval_data[3], &(eigenval_data[3]));
	
	// Calculate Whitening Matrix
	float eigvec_t[4], eigval_inv[4];//, whitening_data[4];
	arm_matrix_instance_f32 eigenvalue_sqrt_matrix = { 2, 2, eigenval_data };
	arm_matrix_instance_f32 eigvec_transpose = { 2, 2, eigvec_t};
	arm_matrix_instance_f32 eigenvec_matrix = {2, 2, eigenvec_data};
	arm_matrix_instance_f32 eigenval_sqrt_inv_matrix = { 2, 2, eigval_inv };
	arm_matrix_instance_f32 whitening_matrix = { 2, 2, whitening_data };
	arm_mat_trans_f32(&eigenvec_matrix, &eigvec_transpose);
	arm_mat_inverse_f32(&eigenvalue_sqrt_matrix, &eigenval_sqrt_inv_matrix);
	arm_mat_mult_f32(&eigenval_sqrt_inv_matrix, &eigvec_transpose, &whitening_matrix);
	
	// Whiten matrix (whitening_matrix * centered_matrix)
	float partial_centered_data[processing_chunk_size / 2];
	float partial_whitened_data[processing_chunk_size / 2];
	arm_matrix_instance_f32 partial_centered_matrix = { 2, (processing_chunk_size / sizeof(float)), partial_centered_data };
	arm_matrix_instance_f32 partial_whitened_matrix = { 2, (processing_chunk_size / sizeof(float)), partial_whitened_data };
	for (index = 0; index < row_size; index += processing_chunk_size) {
		// determine amount of data to read
		size_t data_size;
		if (index + processing_chunk_size < row_size) {
			data_size = processing_chunk_size;
		} else {
			data_size = row_size - index;
		}
		// Get row 0 chunk from FLASH
		BSP_QSPI_Read((uint8_t*)partial_centered_data, pDest + index, data_size);
		// Get row 1 chunk from FLASH
		BSP_QSPI_Read((uint8_t*)&(partial_centered_data[processing_chunk_size/4]), pDest + row_size + index, data_size);
		
		// Compute partial whitened matrix
		arm_mat_mult_f32(&whitening_matrix, &partial_centered_matrix, &partial_whitened_matrix);
		
		// Write back partial whitened matrix
		BSP_QSPI_Write((uint8_t*)partial_whitened_data, pDest2 + index, data_size);
		BSP_QSPI_Write((uint8_t*)&(partial_whitened_data[processing_chunk_size/4]), pDest2 + row_size + index, data_size);
		/*
		if(index == 0){
			for(int i =0;i<10;i++){
				printf("%f  ",partial_whitened_data[i]);
				printf("%f  \n",partial_whitened_data[i+processing_chunk_size/4]);
			}
		}
		
		BSP_QSPI_Read((uint8_t*)partial_whitened_data, pDest2 + index, data_size);
		BSP_QSPI_Read((uint8_t*)&(partial_whitened_data[processing_chunk_size/4]), pDest2 + row_size + index, data_size);
		
		if(index == 0){
			for(int i =0;i<10;i++){
				printf("%f  ",partial_whitened_data[i]);
				printf("%f  \n",partial_whitened_data[i+processing_chunk_size/4]);
			}
		}
		*/
	}
#ifdef STATUS
	printf(" Success\n");
#endif
}

void fast_ica_f32(uint32_t pSrc, uint32_t pDest, size_t size)
{
	float whitening_data[4];
	prewhiten_f32(pSrc, pDest, size,whitening_data);								// prewhiten_f32(pSrc, pDest, size);
	pDest+=size;
	uint32_t samples_per_chunk = processing_chunk_size/sizeof(float);
	const uint32_t row_size = size / 2;
	//fpica
	//added two blocks to store matrix
	uint32_t first_component_location = pDest + size;
	uint32_t second_component_location = first_component_location + row_size;
	//initialize
	int round = 0;
	float wCurrent[2];		//current weight vector
	float wOld[2] = {0,0};	// old weight vector
	float Btrans[4];
	float wCal[2];//temporary matrix used during calculation of w
	float partial_whitened_data[processing_chunk_size / 2];
	arm_matrix_instance_f32 WMat = { 2, 1, wCurrent };
	arm_matrix_instance_f32 BTransMat = { 2, 2, Btrans};
	arm_matrix_instance_f32 wCalMat = {2, 1, wCal };
	
	// Initialize 
	float norm;
	wCurrent[0] = (float)rand();
	wCurrent[1] = (float)rand();
	arm_sqrt_f32(wCurrent[0]*wCurrent[0]+wCurrent[1]*wCurrent[1], &norm);
	wCurrent[0]/=norm;
	wCurrent[1]/=norm;
	
	#ifdef STATUS
	printf("\tCalculating weight vector...");
	#endif
	// Get weight vector to form basis matrix
	for (round = 0; round < MAX_ITERATIONS; round++) {
		// Check if vector is still changing
		//printf("ICA vector = \n%f %f\n\n",wCurrent[0],wCurrent[1]);

		
		wCal[0] = wCurrent[0] - wOld[0];
		wCal[1] = wCurrent[1] - wOld[1];
		arm_sqrt_f32(wCal[0]*wCal[0] + wCal[1]*wCal[1], &norm);
		//printf("norm = %f\n",norm);
		if (norm < EPSILON) break;
		wCal[0] = wCurrent[0] + wOld[0];
		wCal[1] = wCurrent[1] + wOld[1];
		arm_sqrt_f32(wCal[0]*wCal[0] + wCal[1]*wCal[1], &norm);
		//printf("norm = %f\n",norm);
		if (norm < EPSILON) break;
		
		// Update weight
		wOld[0] = wCurrent[0];
		wOld[1] = wCurrent[1];
		
		float inv_num_samples = (float)1 / (row_size / sizeof(float));
		
		// Get weight vector
		float partial_whitened_transpose[processing_chunk_size / 2];
		float partial_vector[processing_chunk_size / 4], temp[processing_chunk_size / 4];
		size_t index;
		float first_partial = 0, second_partial = 0;
		for (index = 0; index < row_size; index += processing_chunk_size) {
			size_t data_size;
			if (index + processing_chunk_size < row_size) {
				data_size = processing_chunk_size;
			} else {
				data_size = row_size - index;
			}
			
			BSP_QSPI_Read((uint8_t*)partial_whitened_data, pDest + index, data_size);
			BSP_QSPI_Read((uint8_t*)&(partial_whitened_data[processing_chunk_size/4]), pDest + row_size + index, data_size);

			arm_matrix_instance_f32 partial_W = { 2, data_size / sizeof(float), partial_whitened_data };
			arm_matrix_instance_f32 partial_W_T = { data_size / sizeof(float), 2, partial_whitened_transpose };
			arm_matrix_instance_f32 partial_vec_mat = { data_size / sizeof(float), 1, partial_vector };
			arm_matrix_instance_f32 temp_mat = { data_size / sizeof(float), 1, temp };
			// W^T * w 
			arm_mat_trans_f32(&partial_W, &partial_W_T);
			arm_mat_mult_f32(&partial_W_T, &WMat, &partial_vec_mat);
			
			/*
			if(index == 0){
				for(int i =0;i<10;i++){
					printf("%f  ",partial_whitened_data[i]);
					printf("%f  ",partial_whitened_data[i+processing_chunk_size/sizeof(float)]);
					printf("%f\n",partial_vector[i]);
				}
			}
			*/
			// (W^T * w)^3 / n
			
			arm_mult_f32(partial_vector, partial_vector, temp, sizeof(partial_vector) / sizeof(float));
			arm_mult_f32(partial_vector, temp, temp, sizeof(partial_vector) / sizeof(float));
			arm_scale_f32(temp, inv_num_samples, temp, sizeof(temp) / sizeof(float));
			/*
			if(index == 0){
				printf("(W^T * w)^3 / n: \n");
				for(int i =0;i<10;i++){
					printf("%f\n",temp[i]);
				}
			}
			*/
			//bug fixed, following parts not tested
			// W * (W^T * w)^3 / n
			arm_mat_mult_f32(&partial_W, &temp_mat, &wCalMat);
			first_partial += wCal[0];
			second_partial += wCal[1];
		}
		// W * (W^T * w)^3 / n - 3w
		wCurrent[0] = first_partial - 3*wCurrent[0];
		wCurrent[1] = second_partial - 3*wCurrent[1];
		arm_sqrt_f32(wCurrent[0]*wCurrent[0]+wCurrent[1]*wCurrent[1], &norm);
		wCurrent[0]/=norm;
		wCurrent[1]/=norm;
	}
	
	// Form transpose of basis matrix
	Btrans[0] = wCurrent[0];
	Btrans[1] = wCurrent[1];
	Btrans[2] = wCurrent[1];
	Btrans[3] = -wCurrent[0];
#ifdef STATUS
	printf(" Success\n\tSeparating components...");
#endif
#ifdef DEBUG	
	printf("whitening matrix = \n%f %f\n%f %f\n\n",whitening_data[0],whitening_data[1],whitening_data[2],whitening_data[3]);
	printf("Btrans matrix = \n%f %f\n%f %f\n\n",Btrans[0],Btrans[1],Btrans[2],Btrans[3]);
#endif
	arm_matrix_instance_f32 whitening_mat = { 2, 2, whitening_data};
	float Wdata[4];
	arm_matrix_instance_f32 W_mat = { 2, 2, Wdata};
	arm_mat_mult_f32(&BTransMat, &whitening_mat, &W_mat);//getting W matrix
#ifdef DEBUG	
	printf("ICA matrix = \n%f %f\n%f %f\n\n",Wdata[0],Wdata[1],Wdata[2],Wdata[3]);
#endif

	
	float partial_component_data[processing_chunk_size / 2];	
	size_t index;
	for (index = 0; index < row_size; index += processing_chunk_size) {
		// Calculate how much data to actually read/write
		size_t data_size;
		if (index + processing_chunk_size < row_size) {
			data_size = processing_chunk_size;
		} else {
			data_size = row_size - index;
		}
		
		arm_matrix_instance_f32 whitened_mat = {2, data_size / sizeof(float), partial_whitened_data};
		arm_matrix_instance_f32 component_mat = {2, data_size / sizeof(float), partial_component_data};
		
		BSP_QSPI_Read((uint8_t*)partial_whitened_data, pSrc + index, data_size);
		BSP_QSPI_Read((uint8_t*)&(partial_whitened_data[processing_chunk_size/4]), pSrc + row_size + index, data_size);
		
		arm_mat_mult_f32(&W_mat, &whitened_mat, &component_mat);
		
		
		BSP_QSPI_Write((uint8_t*)partial_component_data, first_component_location + index, data_size);
		BSP_QSPI_Write((uint8_t*)&(partial_component_data[processing_chunk_size/4]), second_component_location + index, data_size);
		
		
		BSP_QSPI_Read((uint8_t*)partial_component_data, first_component_location + index, data_size);
		BSP_QSPI_Read((uint8_t*)&(partial_component_data[processing_chunk_size/4]), second_component_location + index, data_size);
#ifdef DEBUG	
		if(index == 0){
			for(int i =0;i<30;i++){
				printf("%f  ",partial_whitened_data[i]);
				printf("%f  ",partial_whitened_data[i+processing_chunk_size/sizeof(float)]);
				printf("%f  ",partial_component_data[i]);
				printf("%f\n",partial_component_data[i+processing_chunk_size/sizeof(float)]);
			}
		}
#endif
	}
#ifdef STATUS
	printf(" Success\n");
#endif
}
