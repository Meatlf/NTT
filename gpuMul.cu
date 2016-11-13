/*
 ============================================================================
 Name        : gpuMul.cu
 Author      : jiyang
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA compute reciprocals
 ============================================================================
 */

#include <iostream>
#include <numeric>
#include <stdlib.h>
#include "NTT.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftXt.h>
#include <stdint.h>
#include <stdlib.h>
#include "ModP.h"
#include "kernel.h"
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

using namespace NTL;

void outputGPU(uint64_t* input, uint64_t length){
	uint64_t* temp = (uint64_t*)malloc(sizeof(uint64_t) * length);
	cudaMemcpy(input, temp, length * sizeof(uint64_t), cudaMemcpyDeviceToHost);
	
	for(int i = 0; i < length; i++){
		cout << temp[i] << ", ";
		if(i == length - 1){
			cout << endl;
		}
	}
	free(temp);
	return;
}


void outputGPU(uint32_t* input, uint64_t length){
	uint32_t* temp = (uint32_t*)malloc(sizeof(uint32_t) * length);
	cudaMemcpy(input, temp, length * sizeof(uint32_t), cudaMemcpyDeviceToHost);

	for(int i = 0; i < length; i++){
		cout << temp[i] << ", ";
		if(i == length - 1){
			cout << endl;
		}
	}
	free(temp);
	return;
}

int main(void)
{

	int length = 65536;
	int size = length / 2;

	uint32_t* x = (uint32_t*)malloc(sizeof(uint32_t) * length);
	uint32_t* x_gpu;
	cudaMalloc((void**)&x_gpu, length * sizeof(uint32_t));

	uint32_t* y = (uint32_t*)malloc(sizeof(uint32_t) * length);
	uint32_t* y_gpu;
	cudaMalloc((void**)&y_gpu, length * sizeof(uint32_t));

	uint64_t* X = (uint64_t*)malloc(sizeof(uint32_t) * length);
	uint64_t* X_gpu;
	cudaMalloc((void**)&X_gpu, length * sizeof(uint64_t));


	uint64_t* Y = (uint64_t*)malloc(sizeof(uint64_t) * length);
	uint64_t* Y_gpu;
	cudaMalloc((void**)&Y_gpu, length * sizeof(uint64_t));

	uint32_t* z = (uint32_t*)malloc(sizeof(uint32_t) * length);
	uint64_t* Z_gpu;
	cudaMalloc((void**)&Z_gpu, length * sizeof(uint64_t));
	uint32_t* z_gpu;
	cudaMalloc((void**)&z_gpu, length * sizeof(uint32_t));

	for(int i = 0; i < size; i++){
		x[i] = 16777215;
		y[i] = 16777215;
	}

	cudaMemcpy(x, x_gpu, length * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(y, y_gpu, length * sizeof(uint32_t), cudaMemcpyHostToDevice);

	initNtt(length);

	_ntt(X_gpu, x_gpu, 0, 0, length);

	_ntt(Y_gpu, y_gpu, 0, 0, length);

	dotMul(X_gpu, Y_gpu, Z_gpu, length);

	free(x);free(y);free(X);free(Y);free(z);
	cudaFree(x_gpu);cudaFree(y_gpu);cudaFree(X_gpu);cudaFree(Y_gpu);cudaFree(Z_gpu);cudaFree(z_gpu);

	return 0;
}
