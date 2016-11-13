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
   uint32t tid= threadIdx.x , to , from ;
   to =( tid & 0x07 ) | ( ( tid & 0x38)<<1);
   //2−way bank c o n f l i c t s
   buffer [ to +0x000 ]= samples [ 0 ] ;
   buffer [ to +0x008 ]= samples [ 1 ] ;
   buffer [ to +0x080 ]= samples [ 2 ] ;
   buffer [ to +0x088 ]= samples [ 3 ] ;
   buffer [ to +0x100 ]= samples [ 4 ] ;
   buffer [ to +0x108 ]= samples [ 5 ] ;
   buffer [ to +0x180 ]= samples [ 6 ] ;
   buffer [ to +0x188 ]= samples [ 7 ] ;
   syncthreads ( ) ;
   from =( tid & 0x0F ) | ( ( tid & 0x30)<<3);
   // 2−way bank c o n f l i c t s
   samples [0]= transpose [ from+0x000 ] ;
   samples [1]= transpose [ from+0x010 ] ;
   samples [2]= transpose [ from+0x020 ] ;
   samples [3]= transpose [ from+0x030 ] ;
   samples [4]= transpose [ from+0x040 ] ;
   samples [5]= transpose [ from+0x050 ] ;
   samples [6]= transpose [ from+0x060 ] ;
   samples [7]= transpose [ from+0x070 ] ;
  }
int main(void)
  {
     TransposeSamples<<<1,64>>>();
     return 0;
  }
