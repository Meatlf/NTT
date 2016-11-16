#include <cuda.h>
#include <cuda_runtime.h>
#include <stdint.h>
#include "ModP.h"
#include "kernel.h"

using namespace cuHE;
__global__ void dotMul_kernel (uint64_t* x, uint64_t* y, uint64_t* output, int length){
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	uint64_t _x, _y, _z;

	if(index < length){
		_x = x[index];
		_y = y[index];
		_z = _mul_modP(_x, _y);
		
		output[index] = _z;
	}
}


void dotMul(uint64_t* x, uint64_t*y, uint64_t* output, int length){
	dotMul_kernel<<<length/512 + 1, 512>>>(x, y, output, length);
}
