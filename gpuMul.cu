/*
 ============================================================================
 Name        : gpuMul.cu
 Author      : ttz 
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

		
int main(){
	int len=16;
	uint64 *hx,*dx;

	cudaMallocHost(&hx,len*sizeof(uint64));

	cudaMalloc(&dx,len*sizeof(uint64));

  for(int i=0;i<4;i++){
		hx[i]=11;		
	}
	for(int i=4;i<16;i++){
		hx[i]=33;
	}
  	cudaMemcpy(dx,hx,len*sizeof(uint64),cudaMemcpyHostToDevice);
	ntt_16_1(dx);
	cudaMemcpy(hx,dx,len*sizeof(uint64),cudaMemcpyDeviceToHost);

	for(int i=0;i<len;i++){
		cout<<hx[i]<<endl;
	}
	cudaFree(dx);
	cudaFreeHost(hx);
	return 0;
}
