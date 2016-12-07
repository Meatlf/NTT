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

		
#define len 256

int main(){
	const ZZ P=to_ZZ(0xffffffff00000001);
	const ZZ root_256=to_ZZ((uint64)14041890976876060974);

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
		hx[i]=i;	
	 	ht[i]=0;
	}
  for(int i=0;i<16;i++){
		for(int k=0;k<16;k++){
				conv(h_roots[16*i+k],PowerMod(root_256,i*k,P));
		//test:		cout<<h_roots[16*i+k]<<endl;
		}
	}
  	cudaMemcpy(dx,hx,len*sizeof(uint32),cudaMemcpyHostToDevice);
  	cudaMemcpy(d_roots,h_roots,len*sizeof(uint64),cudaMemcpyHostToDevice);
  	ntt_256(dt,dx,d_roots);
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
