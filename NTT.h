#ifndef _NTT_H_
#define _NTT_H_


#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <cufft.h>
#include "cuda.h"
#include <cuda_runtime.h>
#include <cufftXt.h>
#include <stdint.h>
NTL_CLIENT
float ntt_64K(uint64_t* Y,uint64_t* X,uint32_t* x,uint64_t* root,uint64_t* root_64K);
void ntt_256(uint64_t* X,uint32_t* x,uint64_t* root);
void initNtt(int length) ;

uint64_t **ptrNttSwap();
uint64_t *ptrNttSwap(int dev) ;

void _ntt(uint64_t *X, uint32_t *x, int dev, cudaStream_t st, uint64_t length);
void _intt(uint32_t *x, uint64_t *X, uint32_t crtidx, int dev, cudaStream_t st, uint64_t length);

#endif
