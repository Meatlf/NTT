/*Function：Each transpose operator uses shared memory to pass samples between different warps in the same block.*/
//includes,system
#include<stdio.h>

//includes CUDA Runtime
#include<cuda_runtime.h>

typedef unsigned int uint32;
typedef long int uint64;
__global__ void TransposeSamples()
  {
   __shared__ uint64 buffer[512];
   register uint64 samples[8]={0,1,2,3,4,5,6,7};

   uint32 tid=threadIdx.x;
   uint32 low=tid%8,high=tid/8;
   int index;

   /*write the samples to shared memory*/
   for(index=0;index<8;index++)
	buffer[tid+index*64]=samples[index];
   /*synchronize the threads*/
   __syncthreads();
   /*read the values back using a mapping*/
   for(index=0;index<8;index++)
	samples[index]=buffer[low+index*8+high*64];
  }
int main(void)
  {
     TransposeSamples<<<1,64>>>();
     return 0;
  }

