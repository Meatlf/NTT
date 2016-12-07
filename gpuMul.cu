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
#define len_64K 65536
int main(){
	float time;
	const ZZ P=to_ZZ(0xffffffff00000001);
	const ZZ root_256=to_ZZ((uint64)14041890976876060974);
	const ZZ root_64K=to_ZZ((uint64)15893793146607301539);
	uint32 *hx,*dx;
	uint64 *ht,*dt;
	uint64 *hy,*dy;

	uint64 *h_roots,*d_roots;
  uint64 *h_roots_64K,*d_roots_64K;
  
	cudaMallocHost(&hx,len_64K*sizeof(uint32));
	cudaMallocHost(&ht,len_64K*sizeof(uint64));
	cudaMallocHost(&hy,len_64K*sizeof(uint64));

	cudaMallocHost(&h_roots,len*sizeof(uint64));
  cudaMallocHost(&h_roots_64K,len_64K*sizeof(uint64));

	cudaMalloc(&dx,len_64K*sizeof(uint32));
	cudaMalloc(&dt,len_64K*sizeof(uint64));
	cudaMalloc(&dy,len_64K*sizeof(uint64));

	cudaMalloc(&d_roots,len*sizeof(uint64));
	cudaMalloc(&d_roots_64K,len_64K*sizeof(uint64));
  for(int i=0;i<len_64K;i++){
		hx[i]=i;	
	 	ht[i]=0;
		hy[i]=0;
	}
  for(int i=0;i<16;i++){
		for(int k=0;k<16;k++){
				conv(h_roots[16*i+k],PowerMod(root_256,i*k,P));
		//test:		cout<<h_roots[16*i+k]<<endl;
		}
	}
	for(int i=0;i<256;i++){
		for(int k=0;k<256;k++){
				conv(h_roots_64K[256*i+k],PowerMod(root_64K,i*k,P));
		}
	}
  	cudaMemcpy(dx,hx,len_64K*sizeof(uint32),cudaMemcpyHostToDevice);
  	cudaMemcpy(d_roots,h_roots,len*sizeof(uint64),cudaMemcpyHostToDevice);
		cudaMemcpy(d_roots_64K,h_roots_64K,len_64K*sizeof(uint64),cudaMemcpyHostToDevice);
  	time=ntt_64K(dy,dt,dx,d_roots,d_roots_64K);
  	cudaMemcpy(hy,dy,len_64K*sizeof(uint64),cudaMemcpyDeviceToHost);
	
	cout<<time<<endl;

	for(int i=0;i<len_64K;i++){
		cout<<hy[i]<<endl;
	}
	cudaFree(dx);
	cudaFree(dt);
	cudaFree(dy);
	cudaFree(d_roots);
	cudaFree(d_roots_64K);
	cudaFreeHost(hx);
	cudaFreeHost(ht);
	cudaFreeHost(hy);
	cudaFreeHost(h_roots);
	cudaFreeHost(h_roots_64K);
	return 0;
}
