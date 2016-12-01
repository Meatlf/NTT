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

		
#define len 4096

int main(){
	const ZZ P=to_ZZ(0xffffffff00000001);
	const ZZ root_4K=to_ZZ((uint64)10974926054405199669);

	uint32 *hx,*dx;
	uint64 *ht,*dt;
	uint64 *h_roots,*d_roots;

	cudaMallocHost(&hx,len*sizeof(uint32));
	cudaMallocHost(&ht,len*sizeof(uint64));
	cudaMallocHost(&h_roots,len*sizeof(uint64));

	cudaMalloc(&dx,len*sizeof(uint32));
	cudaMalloc(&dt,len*sizeof(uint64));
	cudaMalloc(&d_roots,len*sizeof(uint64));

  for(int i=0;i<len;i++){
		hx[i]=i;		//32 bits random number
	 	ht[i]=0;
	}
  for(int i=0;i<64;i++){
		for(int k=0;k<64;k++){
				conv(h_roots[64*i+k],PowerMod(root_4K,i*k,P));
		}
	}
  	cudaMemcpy(dx,hx,len*sizeof(uint32),cudaMemcpyHostToDevice);
	cudaMemcpy(d_roots,h_roots,len*sizeof(uint64),cudaMemcpyHostToDevice);
  	ntt_4K_1(dt,dx,d_roots);
	cudaMemcpy(ht,dt,len*sizeof(uint64),cudaMemcpyDeviceToHost);

	for(int i=0;i<len;i++){
		cout<<ht[i]<<endl;
	}
	cudaFree(dx);
	cudaFree(dt);
	cudaFree(d_roots);
	cudaFreeHost(hx);
	cudaFreeHost(ht);
	cudaFreeHost(h_roots);
	return 0;
}
