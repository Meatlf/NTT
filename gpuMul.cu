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
	uint32 *hx,*dx;
	uint64 *ht,*dt;

	cudaMallocHost(&hx,len*sizeof(uint32));
	cudaMallocHost(&ht,len*sizeof(uint64));

	cudaMalloc(&dx,len*sizeof(uint32));
	cudaMalloc(&dt,len*sizeof(uint64));

  for(int i=0;i<len;i++){
		hx[i]=111;		//32 bits random number
	 	ht[i]=0;
		//cout<<"test"<<hx[i];   //for test the function "rand"
	}
  	cudaMemcpy(dx,hx,len*sizeof(uint32),cudaMemcpyHostToDevice);
	ntt_16_1(dt,dx);
	cudaMemcpy(ht,dt,len*sizeof(uint64),cudaMemcpyDeviceToHost);

	for(int i=0;i<len;i++){
		cout<<ht[i]<<endl;
	}
	cudaFree(dx);
	cudaFree(dt);
	cudaFreeHost(hx);
	cudaFreeHost(ht);
	return 0;
}
