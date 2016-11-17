// Copyright (C) IBM, All Rights Reserved

//////#include "StdAfx.h"""
#include <cstring>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <NTL/tools.h>
#include <NTL/ZZ.h>
//#include <xmmintrin.h>
#include <math.h>

NTL_CLIENT

/* Modes of call:
 *   ops = 0: only read key from input (used to test I/O)
 *         1: only generate key
 *         2: generate key and perform re-crypt test
 *
 *   mode = 0: write (read) entire key
 *          1: write (read) w, but not the indexes
 *          2: write (read) only public key, without w or indexes
 ********************************************************************/
#define _N        0xFFFFFFFF00000001ull
#define _ROOT     0xE9653C8DEFA860A9ull    /* 3*2^32th root of unity */
#define _CRITICAL 0xFFFFFFFFul

uint64_t _add(uint64_t a, uint64_t b) {
   uint64_t sum=a+b;

   if(sum<a) {
	if(sum>=_N)
	 sum=sum-_N-_N;
	else
	 sum=sum-_N;
   }
   else if(sum>=a && sum>=_N)
	sum=sum-_N;
   return sum;
}

uint64_t _subtract(uint64_t a, uint64_t b) {
   uint64_t sum;

   if(b>=_N)
    b-=_N;
   if(a>=b)
    sum=a-b;
   else
    sum=a-b+_N;
   return sum;
}

uint64_t _multiply(uint64_t a, uint64_t b) {
   uint64_t al=(uint32_t)a, bl=(uint32_t)b, ah=a>>32, bh=b>>32;
   uint64_t albl=al*bl, albh=al*bh, ahbl=ah*bl, ahbh=ah*bh;
   uint64_t upper, carry=0;
   uint32_t uu, ul;

   upper=(albl>>32)+albh+ahbl;
   if(upper<ahbl)
    carry=0x100000000ull;
   upper=(upper>>32)+ahbh+carry;

   uu=upper>>32;
   ul=upper;

   if(ul==0)
    upper=_N-uu;
   else
    upper=(upper<<32)-uu-ul;
   //printf("upper=%Lx\n", upper>>32);
   //printf("a*b=%Lx\n", a*b);

   return _add(a*b, upper);
}

uint64_t _power(uint64_t a, uint64_t k) {
   uint64_t current=1, square=a;

   while(k>0) {
	if((k&1)!=0)
	 current=_multiply(current, square);
	square=_multiply(square, square);
	k=k>>1;
   }
   return current;
}

uint64_t _normalize(uint64_t a) {
   if(a>=_N)
    return a-_N;
   return a;
}

void extendedGCD(uint64_t a, uint64_t b, int64_t *s, int64_t *t) {
   int64_t x, y;

   if(a%b==0) {
    *s=0;
    *t=1;
   }
   else {
    extendedGCD(b, a%b, &x, &y);
    *s=y;
    *t=x-y*(a/b);
   }
}

uint64_t _inverse(uint64_t a) {
   int64_t s, t;

   extendedGCD(_N, a, &s, &t);
   if(t<0)
    return t+_N;
   return t;
}

uint64_t _root(uint64_t size) {
   uint64_t k=(3ull<<32)/size;

   return _normalize(_power(_ROOT, k));
}

uint64_t _inverseRoot(uint64_t size) {
   uint64_t k=(3ull<<32)/size;

   return _normalize(_power(_ROOT, (3ull<<32)-k));
}

void transpose(uint64_t *x, uint64_t *X, uint32_t size, uint32_t xLength) {
   int group, index;

   for(group=0;group<size/xLength;group++)
	for(index=0;index<xLength;index++)
     X[group*xLength+index]=x[index*size/xLength+group];
}

void smallFFTNoMultiply(uint64_t *x, uint64_t *X, uint32_t size, int inverse) {
   int      i, j;
   uint64_t total, r, value;

   if(inverse)
    r=_inverseRoot(size);
   else
    r=_root(size);
   for(i=0;i<size;i++) {
	total=0;
	for(j=0;j<size;j++) {
     value=_multiply(_power(r, i*j), x[j]);
	 total=_add(total, value);
    }
	X[i]=_normalize(total);
   }
}

void smallFFT(uint64_t *x, uint64_t *X, uint32_t size, int inverse) {
   int      i, j;
   uint64_t total, r, inv=1, value;

   if(inverse) {
    r=_inverseRoot(size);
    inv=_inverse(size);
   }
   else
    r=_root(size);
   //printf("root8=%d\n", r);
   for(i=0;i<size;i++) {
	total=0;
	for(j=0;j<size;j++) {
     value=_multiply(_power(r, i*j), x[j]);
	 total=_add(total, value);
    }
	X[i]=_normalize(_multiply(total, inv));
   }
   //printf("smallFFT:X[0]=%x\n",X[0]);
}

void split_by2(uint32_t *x, uint32_t xLength, uint64_t *X, uint32_t XLength) {
   int      index, byteCount=0, current=0;
   uint32_t value=0;

   for(index=0;index<XLength;index++) {
	  if(byteCount==0 && current<xLength) {
	     value=x[current++];
	     byteCount=4;
      }
      X[index]=value & 0xFFFF;
      byteCount-=2;
      value=value>>16;
   }
}

void split_by3(uint32_t *x, uint32_t xLength, uint64_t *X, uint32_t XLength) {
   int      byteCount=0, current=0, index;
   uint64_t value=0;

   for(index=0;index<XLength;index++) {
    if(byteCount<3) {
     if(current<xLength)
      value+=((uint64_t)x[current++])<<byteCount*8;
	 byteCount+=4;
    }
	X[index]=value & 0xFFFFFF;
	byteCount-=3;
	value=value>>24;
   }
}

void resolve_by2(uint64_t *x, uint32_t length) {
   uint64_t carry=0, current;
   int      index;

   for(index=0;index<length;index++) {
	  current=x[index];
	  carry+=(current & 0xFFFF);
	  x[index]=carry & 0xFFFF;
	  carry=(carry>>16) + (current>>16);
   }
   if(carry>0) {
	  printf("resolve_by2 failed! abort!\n");
	  //exit(1);
   }
}

void resolve_by3(uint64_t *x, uint32_t length) {
  uint64_t carry = 0, current;
  int index;

  for(index = 0;index < length; index++) {
    current = x[index];
    carry += (current & 0xFFFFFF);
    x[index] = carry & 0xFFFFFF;
    carry = (carry >> 24) + (current >> 24);
  }

  if(carry>0) {
    printf("resolve_by3 failed! abort!\n");
    exit(1);
  }
}

void join_by2(uint64_t *x, uint32_t xLength, uint32_t *X, uint32_t XLength) {
   int byteCount = 0, current = 0, index;
   uint64_t value = 0;

  value = 0;
  byteCount = 0;
  current = 0;
  for(index = 0;index < XLength; index++) {
    while(byteCount<4) {
      if(current<xLength)
        value+=(x[current++])<<byteCount*8;
      byteCount+=2;
    }
      X[index]=value;
      value=value>>32;
      byteCount-=4;
   }
}

void join_by3(uint64_t *x, uint32_t xLength, uint32_t *X, uint32_t XLength) {
   int      byteCount=0, current=0, index;
   uint64_t value=0;

   value=0;
   byteCount=0;
   current=0;
   for(index=0;index<XLength;index++) {
    while(byteCount<4) {
     if(current<xLength)
      value+=(x[current++])<<byteCount*8;
     byteCount+=3;
	}
    X[index]=value;
    value=value>>32;
    byteCount-=4;
   }
}

void fft64AndTwiddle(uint64_t *x, uint64_t *X, uint64_t r, int inverse) {
   uint64_t smallX[64], smallY[64];
   int      fft, index;

   for(fft=0;fft<8;fft++) {
	for(index=0;index<64;index++)
	 smallX[index]=x[index*8+fft];
	smallFFT(smallX, smallY, 64, inverse);
	for(index=0;index<64;index++) {
     if(inverse)
      X[index*8+fft]=_multiply(_multiply(smallY[index], _power(r, index)), 64);
     else
      X[index*8+fft]=_multiply(smallY[index], _power(r, index));
    }
   }
}

void load64(uint64_t *from, uint64_t *x, uint32_t offset, uint32_t stride) {
   uint32_t group, index;

   for(group=0;group<64;group++) {
    for(index=0;index<8;index++)
     x[group*8+index]=from[offset+group*stride+index];
   }
}

void store64(uint64_t *to, uint64_t *x, uint32_t offset, uint32_t stride) {
   uint32_t group, index;

   for(group=0;group<64;group++) {
    for(index=0;index<8;index++)
     to[offset+group*stride+index]=x[group*8+index];
   }
}

void largeFFT(uint64_t *x, uint64_t *X, uint32_t size, int inverse) {
   int      count=1, fft, i, j, k, l, m;
   uint64_t r, inv=1, omega, currentOmega, c0, c1;
   uint64_t *buffer, *to, *from, *swap, *smallX, *smallY;

   if(size%3==0) {
		size=size/3;
		count=3;
   }
   if(inverse) {
		r=_inverseRoot(size);
		inv=_inverse(size);
   }
   else
		r=_root(size); //
		
   buffer=(uint64_t *)malloc(sizeof(uint64_t)*size);

   for(fft=0;fft<count;fft++) {
		for(i=0;i<size;i++)
			buffer[i]=x[count*i+fft];
		from=buffer;
		to=X+fft*size;
		l=size/2;
		m=1;
		while(l>=1) {
			omega=_power(r, size/(2*l));
			currentOmega=1;
			for(j=0;j<l;j++) {
				for(k=0;k<m;k++) {
					c0=from[k+j*m];
					c1=from[k+j*m+l*m];
					to[k+2*j*m]=_add(c0, c1);
					to[k+2*j*m+m]=_multiply(currentOmega, _subtract(c0, c1));
				}
				currentOmega=_multiply(currentOmega, omega);
			}
			swap=from;
			from=to;
			to=swap;
			l=l>>1;
			m=m<<1;
		}
		if(from!=X+fft*size) {
			for(i=0;i<size;i++)
				X[fft*size+i]=from[i];
		}
	}

	if(count>1) {
		r=_root(size*count);
		smallX=(uint64_t *)malloc(sizeof(uint64_t)*count);
		smallY=(uint64_t *)malloc(sizeof(uint64_t)*count);
		for(fft=0;fft<size;fft++) {
			smallX[0]=X[fft];
			smallX[1]=_multiply(_power(r, fft), X[fft+size]);
			smallX[2]=_multiply(_power(r, fft+fft), X[fft+size+size]);
			smallFFT(smallX, smallY, 3, inverse);
			X[fft]=smallY[0];
			X[fft+size]=smallY[1];
			X[fft+size+size]=smallY[2];
		}
		free(smallX);
		free(smallY);
	}

	for(i=0;i<size*count;i++)
		X[i]=_normalize(_multiply(X[i], inv));
	free(buffer);
}

void radix64V1(uint64_t *from, uint64_t *to, uint32_t size, uint32_t l, int inverse) {
   uint64_t smallX[64*8], smallY[64*8];
   uint32_t blockIdxX, gridDimX=64;
   uint32_t perBlock, z, m, j, k, base, current, offset;
   uint32_t fromOffset, fromStride, toOffset, toStride;
   uint64_t r;

   // size=xy=64^n*8*z   x=64^n, y=8*z   1<=z<64 or z=96

   z=size/8;
   while((z%64)==0)
    z=z/64;

   perBlock=size/512/gridDimX;
   m=size/512/z/l;

   for(blockIdxX=0;blockIdxX<gridDimX;blockIdxX++) {
    base=blockIdxX*perBlock;
    current=perBlock;
    while(current>0) {
     current--;

     offset=(base+current)%z;
     k=(base+current)/z/l;
     j=(base+current)/z%l;

     fromOffset=(j+64*k*l)*z*8 + offset*8;
     fromStride=l*z*8;
     toOffset=(j+k*l)*z*8 + offset*8;
     toStride=l*m*z*8;

     if(inverse)
      r=_power(_inverseRoot(64*l), j);
     else
      r=_power(_root(64*l), j);
     load64(from, smallX, fromOffset, fromStride);
     fft64AndTwiddle(smallX, smallY, r, inverse);
     store64(to, smallY, toOffset, toStride);
    }
   }
}

void radix64V2(uint64_t *from, uint64_t *to, uint32_t size, uint32_t l, int inverse) {
   uint64_t smallX[64*8], smallY[64*8];
   uint32_t blockIdxX, gridDimX=64;
   uint32_t perBlock, z, m, j, k, base, current, offset;
   uint32_t fromOffset, fromStride, toOffset, toStride;
   uint64_t r;

   // size=xy=64^n*8*z   x=64^n, y=8*z   1<=z<64 or z=96

   z=size/8;
   while((z%64)==0)
    z=z/64;

   perBlock=size/512/gridDimX;
   m=size/512/z/l;

   printf("m=%d\n",m);

   for(blockIdxX=0;blockIdxX<gridDimX;blockIdxX++) {
    base=blockIdxX*perBlock;
    current=perBlock;
    while(current>0) {
     current--;

     offset=(base+current)%z;
     j=(base+current)/z/m;
     k=(base+current)/z%m;

     fromOffset=(k+j*m)*z*8 + offset*8;
     fromStride=l*m*z*8;
     toOffset=(k+64*j*m)*z*8 + offset*8;
     toStride=m*z*8;

     if(inverse)
      r=_power(_inverseRoot(64*l), j);
     else
      r=_power(_root(64*l), j);
     load64(from, smallX, fromOffset, fromStride);
     fft64AndTwiddle(smallX, smallY, r, inverse);
     store64(to, smallY, toOffset, toStride);
    }
   }
}

void radix64FFT(uint64_t *x, uint64_t *X, uint32_t size, uint32_t level, int inverse) {
   
   uint32_t l, m, j, k, base, current, offset;
   uint64_t inv=1, omega, currentOmega, c0, c1;
   uint64_t *buffer, *to, *from, *swap;
   uint64_t r;
   uint64_t smallX[64], smallY[64];
   uint64_t temp0=0;
   
   buffer=(uint64_t *)malloc(sizeof(uint64_t)*size);

   for(int t=0;t<size;t++)
		buffer[t]=x[t];
   from=buffer;
   to=X;
   l=size/64;
   m=1;

   for(int i=0; i<level; i++) {
	    for(j=0;j<l;j++) {
            if(inverse)
                r=_power(_inverseRoot(64*l), j);
            else
                r=_power(_root(64*l), j);
			for(k=0;k<m;k++) {
				for(int t=0; t<64; t++)
					smallX[t]=from[k+j*m+t*l*m];
				smallFFT(smallX, smallY, 64, inverse);
				for(int t=0; t<64; t++) {
					if(inverse)
						to[k+64*j*m+t*m]=_multiply(smallY[t], _power(r, t));
					else
						to[k+64*j*m+t*m]=_multiply(smallY[t], _power(r, t));
				}
			} 
		}
		swap=from;
		from=to;
		to=swap;
		l=l/64;
		m=m*64;
   }

   for(int i=0;i<size;i++)
       X[i]=from[i];
}

void smallFFT16(uint64_t *x, uint64_t *X, uint32_t size, int inverse) {
   int      i, j;
   uint64_t total, r, inv=1, value;

   if(inverse) {
    r=_inverseRoot(size);
    inv=_inverse(size);
   }
   else
    r=_root(size);
   //printf("root8=%d\n", r);
   for(i=0;i<size;i++) {
	total=0;
	for(j=0;j<size;j++) {
     value=_multiply(_power(r, i*j), x[j]);
	 total=_add(total, value);
    }
	X[i]=_normalize(_multiply(total, inv));
   }
   printf("smallFFT:inv=%016llX\n",inv);
}

void radix16FFT(uint64_t *x, uint64_t *X, uint32_t size, uint32_t level, int inverse) {
   
   uint32_t l, m, j, k, base, current, offset;
   uint64_t inv=1, omega, currentOmega, c0, c1;
   uint64_t *buffer, *to, *from, *swap;
   uint64_t r;
   uint64_t smallX[16], smallY[16];
   uint64_t temp0=0;
   
   buffer=(uint64_t *)malloc(sizeof(uint64_t)*size);

   for(int t=0;t<size;t++)
		buffer[t]=x[t];
   from=buffer;
   to=X;
   l=size/16;
   m=1;
   
   for(int i=0; i<level; i++) {
	 FILE *fp_dataOut = NULL;
	 fp_dataOut = fopen("dataStageOut.txt", "wt");
	    for(j=0;j<l;j++) {
            if(inverse)
                r=_power(_inverseRoot(16*l), j);
            else
                r=_power(_root(16*l), j);
			for(k=0;k<m;k++) {
				for(int t=0; t<16; t++)
					smallX[t]=from[k+j*m+t*l*m];
				//for(int t=0; t<16; t++) {
				//	printf("index=%d, X = %016llX\n", k+j*m+t*l*m, smallX[t]);
				//}
				smallFFT16(smallX, smallY, 16, inverse);
				for(int t=0; t<16; t++) {
					if(inverse)
						to[k+16*j*m+t*m]=_multiply(smallY[t], _power(r, t));
					else
						to[k+16*j*m+t*m]=_multiply(smallY[t], _power(r, t)); //
				}
				for(int t=0; t<16; t++) {
					//printf("index=%d, smallY=%016llX Pow=%016llX Y = %016llX\n", k+j*m+t*l*m, smallY[t], _power(r, t), _multiply(smallY[t], _power(r, t)));
					fprintf(fp_dataOut, "index=%d, smallX=%016llX smallY=%016llX Pow=%016llX Y = %016llX\n", k+j*m+t*l*m,smallX[t], smallY[t], _power(r, t), _multiply(smallY[t], _power(r, t)));
				}
			} 
		}


   //for(int i=0; i<size; i++) {
   //    fprintf(fp_dataOut,"%016llX\n", to[i]);
   //}
   fclose(fp_dataOut);
		swap=from;
		from=to;
		to=swap;
		l=l/16;
		m=m*16;
   }





   for(int i=0;i<size;i++)
       X[i]=from[i];
}

int format_binary(uint64_t x, uint32_t *s)
{
    for (int z = 0; z < 64; z++) {
        s[63-z] = ((x>>z) & 0x0000000000000001) ? 1 : 0;
    }

    return 0;
}


void main()
{
	
	uint32_t  words, samples, index;
	uint32_t  *a, *b, *c;
	uint64_t  *a64, *b64;
	uint64_t  *aFFT, *bFFT, *cFFT, *cFFT1;
	uint64_t  *a64_s0, *a64_s1, *a64_s2, *a64_s3, *a64_s4, *a64_s5;
	uint64_t  *a64_s6, *a64_s7, *a64_s8, *a64_s9, *a64_s10, *a64_s11;
	uint64_t  *a64_s12, *a64_s13, *a64_s14, *a64_s15;

	uint64_t *rom0, *rom1, *rom2, *rom3, *rom4, *rom5, *rom6, *rom7, *rom8;
	uint64_t *rom9, *rom10, *rom11, *rom12, *rom13, *rom14, *rom15;

	uint32_t   d3, d2, d1, d0;

	ZZ        aZZ, bZZ, cZZ;
	uint64_t  r;

	uint64_t inv_64K;
	inv_64K=_inverse(4096);
	printf("inv_64K = %016llX\n", inv_64K);

	//r=_root(64);
	//printf("root64=%d\n",r);
	//r=_inverseRoot(128);
	//printf("root128=%d\n",r);
	//r=_inverseRoot(256);
	//printf("root256=%d\n",r);

	//samples = 24576; //24576*32bits = 786Kbits 
	words = 4096/2*24/32; //samples/4;
	int FFT_SIZE = 4096;//4*words; 64K FFT
	//int FFT_SIZE = 16;
	//words = FFT_SIZE/4;
	int dataIndex;
	int bankIndex;
	
	a=(uint32_t *)malloc(sizeof(uint32_t)*words);
	b=(uint32_t *)malloc(sizeof(uint32_t)*words);
	c=(uint32_t *)malloc(sizeof(uint32_t)*words*2);
	a64=(uint64_t *)malloc(sizeof(uint64_t)*FFT_SIZE);
	b64=(uint64_t *)malloc(sizeof(uint64_t)*FFT_SIZE);
	aFFT=(uint64_t *)malloc(sizeof(uint64_t)*FFT_SIZE);
	bFFT=(uint64_t *)malloc(sizeof(uint64_t)*FFT_SIZE);
	cFFT=(uint64_t *)malloc(sizeof(uint64_t)*FFT_SIZE);
	cFFT1=(uint64_t *)malloc(sizeof(uint64_t)*FFT_SIZE);

	a64_s0 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s1 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s2 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s3 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s4 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s5 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s6 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s7 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s8 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s9 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s10 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s11 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s12 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s13 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s14 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	a64_s15 = (uint64_t *)malloc(sizeof(uint64_t)*256);

	rom0 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom1 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom2 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom3 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom4 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom5 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom6 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom7 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom8 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom9 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom10 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom11 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom12 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom13 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom14 = (uint64_t *)malloc(sizeof(uint64_t)*256);
	rom15 = (uint64_t *)malloc(sizeof(uint64_t)*256);

    for(index=0;index<words;index++) {  //words
	    a[index]=(rand()<<16)^rand();
    }

    for(index=0;index<words;index++) {
	    b[index]=(rand()<<16)^rand();
    }

	split_by3(a, words, a64, FFT_SIZE);
	split_by3(b, words, b64, FFT_SIZE);

	for(int i=0; i<FFT_SIZE; i++) {
		d0 = i % 16;
		d1 = (i>>4) % 16;
		d2 = (i>>8) % 16;
		d3 = (i>>12) % 16;

		dataIndex = i >> 4;
		bankIndex = (d3+d2+d1+d0) % 16;
		switch(bankIndex){
			case 0:a64_s0[dataIndex] = a64[i];break; //a64[i]
			case 1:a64_s1[dataIndex] = a64[i];break; //a64[i]
			case 2:a64_s2[dataIndex] = a64[i];break; //a64[i]
			case 3:a64_s3[dataIndex] = a64[i];break; //a64[i]
			case 4:a64_s4[dataIndex] = a64[i];break; //a64[i]
			case 5:a64_s5[dataIndex] = a64[i];break; //a64[i]
			case 6:a64_s6[dataIndex] = a64[i];break; //a64[i]
			case 7:a64_s7[dataIndex] = a64[i];break; //a64[i]
			case 8:a64_s8[dataIndex] = a64[i];break; //a64[i]
			case 9:a64_s9[dataIndex] = a64[i];break;
			case 10:a64_s10[dataIndex] = a64[i];break;
			case 11:a64_s11[dataIndex] = a64[i];break;
			case 12:a64_s12[dataIndex] = a64[i];break;
			case 13:a64_s13[dataIndex] = a64[i];break;
			case 14:a64_s14[dataIndex] = a64[i];break;
			case 15:a64_s15[dataIndex] = a64[i];break;
		}
	}

	FILE *fp_data = NULL;
	fp_data = fopen("dataIn.txt", "wt");
	for(int i=0; i<FFT_SIZE; i++) {
		fprintf(fp_data,"%x\n", a64[i]);
	}
	fclose(fp_data);

	FILE *fp_data0 = NULL;
	FILE *fp_data1 = NULL;
	FILE *fp_data2 = NULL;
	FILE *fp_data3 = NULL;
	FILE *fp_data4 = NULL;
	FILE *fp_data5 = NULL;
	FILE *fp_data6 = NULL;
	FILE *fp_data7 = NULL;
	FILE *fp_data8 = NULL;
	FILE *fp_data9 = NULL;
	FILE *fp_data10 = NULL;
	FILE *fp_data11 = NULL;
	FILE *fp_data12 = NULL;
	FILE *fp_data13 = NULL;
	FILE *fp_data14 = NULL;
	FILE *fp_data15 = NULL;
	fp_data0 = fopen("dataIn_0.txt", "wt");
	fp_data1 = fopen("dataIn_1.txt", "wt");
	fp_data2 = fopen("dataIn_2.txt", "wt");
	fp_data3 = fopen("dataIn_3.txt", "wt");
	fp_data4 = fopen("dataIn_4.txt", "wt");
	fp_data5 = fopen("dataIn_5.txt", "wt");
	fp_data6 = fopen("dataIn_6.txt", "wt");
	fp_data7 = fopen("dataIn_7.txt", "wt");
	fp_data8 = fopen("dataIn_8.txt", "wt");
	fp_data9 = fopen("dataIn_9.txt", "wt");
	fp_data10 = fopen("dataIn_10.txt", "wt");
	fp_data11 = fopen("dataIn_11.txt", "wt");
	fp_data12 = fopen("dataIn_12.txt", "wt");
	fp_data13 = fopen("dataIn_13.txt", "wt");
	fp_data14 = fopen("dataIn_14.txt", "wt");
	fp_data15 = fopen("dataIn_15.txt", "wt");
	for(int i=0; i<256; i++) {
		fprintf(fp_data0,"%016llX\n", a64_s0[i]);
		fprintf(fp_data1,"%016llX\n", a64_s1[i]);
		fprintf(fp_data2,"%016llX\n", a64_s2[i]);
		fprintf(fp_data3,"%016llX\n", a64_s3[i]);
		fprintf(fp_data4,"%016llX\n", a64_s4[i]);
		fprintf(fp_data5,"%016llX\n", a64_s5[i]);
		fprintf(fp_data6,"%016llX\n", a64_s6[i]);
		fprintf(fp_data7,"%016llX\n", a64_s7[i]);
		fprintf(fp_data8,"%016llX\n", a64_s8[i]);
		fprintf(fp_data9,"%016llX\n", a64_s9[i]);
		fprintf(fp_data10,"%016llX\n", a64_s10[i]);
		fprintf(fp_data11,"%016llX\n", a64_s11[i]);
		fprintf(fp_data12,"%016llX\n", a64_s12[i]);
		fprintf(fp_data13,"%016llX\n", a64_s13[i]);
		fprintf(fp_data14,"%016llX\n", a64_s14[i]);
		fprintf(fp_data15,"%016llX\n", a64_s15[i]);
	}
	fclose(fp_data0);
	fclose(fp_data1);
	fclose(fp_data2);
	fclose(fp_data3);
	fclose(fp_data4);
	fclose(fp_data5);
	fclose(fp_data6);
	fclose(fp_data7);
	fclose(fp_data8);
	fclose(fp_data9);
	fclose(fp_data10);
	fclose(fp_data11);
	fclose(fp_data12);
	fclose(fp_data13);
	fclose(fp_data14);
	fclose(fp_data15);
	
		for(int i=0; i<FFT_SIZE; i++) {
		d0 = i % 16;
		d1 = (i>>4) % 16;
		d2 = (i>>8) % 16;
		d3 = (i>>12) % 16;

		dataIndex = i >> 4;
		bankIndex = (d3+d2+d1+d0) % 16;
		switch(bankIndex){
			case 0:a64_s0[dataIndex] = b64[i];break; //a64[i]
			case 1:a64_s1[dataIndex] = b64[i];break; //a64[i]
			case 2:a64_s2[dataIndex] = b64[i];break; //a64[i]
			case 3:a64_s3[dataIndex] = b64[i];break; //a64[i]
			case 4:a64_s4[dataIndex] = b64[i];break; //a64[i]
			case 5:a64_s5[dataIndex] = b64[i];break; //a64[i]
			case 6:a64_s6[dataIndex] = b64[i];break; //a64[i]
			case 7:a64_s7[dataIndex] = b64[i];break; //a64[i]
			case 8:a64_s8[dataIndex] = b64[i];break; //a64[i]
			case 9:a64_s9[dataIndex] = b64[i];break;
			case 10:a64_s10[dataIndex] = b64[i];break;
			case 11:a64_s11[dataIndex] = b64[i];break;
			case 12:a64_s12[dataIndex] = b64[i];break;
			case 13:a64_s13[dataIndex] = b64[i];break;
			case 14:a64_s14[dataIndex] = b64[i];break;
			case 15:a64_s15[dataIndex] = b64[i];break;
		}
	}

	fp_data = fopen("dataInB.txt", "wt");
	for(int i=0; i<FFT_SIZE; i++) {
		fprintf(fp_data,"%x\n", a64[i]);
	}
	fclose(fp_data);

	fp_data0 = fopen("dataInB_0.txt", "wt");
	fp_data1 = fopen("dataInB_1.txt", "wt");
	fp_data2 = fopen("dataInB_2.txt", "wt");
	fp_data3 = fopen("dataInB_3.txt", "wt");
	fp_data4 = fopen("dataInB_4.txt", "wt");
	fp_data5 = fopen("dataInB_5.txt", "wt");
	fp_data6 = fopen("dataInB_6.txt", "wt");
	fp_data7 = fopen("dataInB_7.txt", "wt");
	fp_data8 = fopen("dataInB_8.txt", "wt");
	fp_data9 = fopen("dataInB_9.txt", "wt");
	fp_data10 = fopen("dataInB_10.txt", "wt");
	fp_data11 = fopen("dataInB_11.txt", "wt");
	fp_data12 = fopen("dataInB_12.txt", "wt");
	fp_data13 = fopen("dataInB_13.txt", "wt");
	fp_data14 = fopen("dataInB_14.txt", "wt");
	fp_data15 = fopen("dataInB_15.txt", "wt");
	for(int i=0; i<256; i++) {
		fprintf(fp_data0,"%016llX\n", a64_s0[i]);
		fprintf(fp_data1,"%016llX\n", a64_s1[i]);
		fprintf(fp_data2,"%016llX\n", a64_s2[i]);
		fprintf(fp_data3,"%016llX\n", a64_s3[i]);
		fprintf(fp_data4,"%016llX\n", a64_s4[i]);
		fprintf(fp_data5,"%016llX\n", a64_s5[i]);
		fprintf(fp_data6,"%016llX\n", a64_s6[i]);
		fprintf(fp_data7,"%016llX\n", a64_s7[i]);
		fprintf(fp_data8,"%016llX\n", a64_s8[i]);
		fprintf(fp_data9,"%016llX\n", a64_s9[i]);
		fprintf(fp_data10,"%016llX\n", a64_s10[i]);
		fprintf(fp_data11,"%016llX\n", a64_s11[i]);
		fprintf(fp_data12,"%016llX\n", a64_s12[i]);
		fprintf(fp_data13,"%016llX\n", a64_s13[i]);
		fprintf(fp_data14,"%016llX\n", a64_s14[i]);
		fprintf(fp_data15,"%016llX\n", a64_s15[i]);
	}
	fclose(fp_data0);
	fclose(fp_data1);
	fclose(fp_data2);
	fclose(fp_data3);
	fclose(fp_data4);
	fclose(fp_data5);
	fclose(fp_data6);
	fclose(fp_data7);
	fclose(fp_data8);
	fclose(fp_data9);
	fclose(fp_data10);
	fclose(fp_data11);
	fclose(fp_data12);
	fclose(fp_data13);
	fclose(fp_data14);
	fclose(fp_data15);

	uint32_t s[16];
	uint64_t root64K=_root(FFT_SIZE);
	FILE *fp_rom1 = NULL;
	FILE *fp_rom2 = NULL;
	FILE *fp_rom3 = NULL;
	FILE *fp_rom4 = NULL;
	FILE *fp_rom5 = NULL;
	FILE *fp_rom6 = NULL;
	FILE *fp_rom7 = NULL;
	FILE *fp_rom8 = NULL;
	FILE *fp_rom9 = NULL;
	FILE *fp_rom10 = NULL;
	FILE *fp_rom11 = NULL;
	FILE *fp_rom12 = NULL;
	FILE *fp_rom13 = NULL;
	FILE *fp_rom14 = NULL;
	FILE *fp_rom15 = NULL;
	fp_rom1 = fopen("ROM1.rcf", "wt");
	fp_rom2 = fopen("ROM2.rcf", "wt");
	fp_rom3 = fopen("ROM3.rcf", "wt");
	fp_rom4 = fopen("ROM4.rcf", "wt");
	fp_rom5 = fopen("ROM5.rcf", "wt");
	fp_rom6 = fopen("ROM6.rcf", "wt");
	fp_rom7 = fopen("ROM7.rcf", "wt");
	fp_rom8 = fopen("ROM8.rcf", "wt");
	fp_rom9 = fopen("ROM9.rcf", "wt");
	fp_rom10 = fopen("ROM10.rcf", "wt");
	fp_rom11 = fopen("ROM11.rcf", "wt");
	fp_rom12 = fopen("ROM12.rcf", "wt");
	fp_rom13 = fopen("ROM13.rcf", "wt");
	fp_rom14 = fopen("ROM14.rcf", "wt");
	fp_rom15 = fopen("ROM15.rcf", "wt");
	for(int i=0; i<256; i++) {
		/*
		fprintf(fp_rom1,"%016llX\n", _power(root64K, i));
		fprintf(fp_rom2,"%016llX\n", _power(root64K, i*2));
		fprintf(fp_rom3,"%016llX\n", _power(root64K, i*3));
		fprintf(fp_rom4,"%016llX\n", _power(root64K, i*4));
		fprintf(fp_rom5,"%016llX\n", _power(root64K, i*5));
		fprintf(fp_rom6,"%016llX\n", _power(root64K, i*6));
		fprintf(fp_rom7,"%016llX\n", _power(root64K, i*7));
		fprintf(fp_rom8,"%016llX\n", _power(root64K, i*8));
		fprintf(fp_rom9,"%016llX\n", _power(root64K, i*9));
		fprintf(fp_rom10,"%016llX\n", _power(root64K, i*10));
		fprintf(fp_rom11,"%016llX\n", _power(root64K, i*11));
		fprintf(fp_rom12,"%016llX\n", _power(root64K, i*12));
		fprintf(fp_rom13,"%016llX\n", _power(root64K, i*13));
		fprintf(fp_rom14,"%016llX\n", _power(root64K, i*14));
		fprintf(fp_rom15,"%016llX\n", _power(root64K, i*15));
		*/
		format_binary(_power(root64K, i), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom1,"%d", s[j]);
		}
		fprintf(fp_rom1,"%d\n", s[15]);

		format_binary(_power(root64K, i*2), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom2,"%d", s[j]);
		}
		fprintf(fp_rom2,"%d\n", s[15]);

		format_binary(_power(root64K, i*3), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom3,"%d", s[j]);
		}
		fprintf(fp_rom3,"%d\n", s[15]);

		format_binary(_power(root64K, i*4), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom4,"%d", s[j]);
		}
		fprintf(fp_rom4,"%d\n", s[15]);

		format_binary(_power(root64K, i*5), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom5,"%d", s[j]);
		}
		fprintf(fp_rom5,"%d\n", s[15]);

		format_binary(_power(root64K, i*6), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom6,"%d", s[j]);
		}
		fprintf(fp_rom6,"%d\n", s[15]);

		format_binary(_power(root64K, i*7), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom7,"%d", s[j]);
		}
		fprintf(fp_rom7,"%d\n", s[15]);

		format_binary(_power(root64K, i*8), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom8,"%d", s[j]);
		}
		fprintf(fp_rom8,"%d\n", s[15]);

		format_binary(_power(root64K, i*9), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom9,"%d", s[j]);
		}
		fprintf(fp_rom9,"%d\n", s[15]);

		format_binary(_power(root64K, i*10), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom10,"%d", s[j]);
		}
		fprintf(fp_rom10,"%d\n", s[15]);

		format_binary(_power(root64K, i*11), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom11,"%d", s[j]);
		}
		fprintf(fp_rom11,"%d\n", s[15]);

		format_binary(_power(root64K, i*12), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom12,"%d", s[j]);
		}
		fprintf(fp_rom12,"%d\n", s[15]);

		format_binary(_power(root64K, i*13), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom13,"%d", s[j]);
		}
		fprintf(fp_rom13,"%d\n", s[15]);

		format_binary(_power(root64K, i*14), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom14,"%d", s[j]);
		}
		fprintf(fp_rom14,"%d\n", s[15]);

		format_binary(_power(root64K, i*15), s);
		for(int j=0; j<15; j++) {
			fprintf(fp_rom15,"%d", s[j]);
		}
		fprintf(fp_rom15,"%d\n", s[15]);


	}
	fclose(fp_rom1);
	fclose(fp_rom2);
	fclose(fp_rom3);
	fclose(fp_rom4);
	fclose(fp_rom5);
	fclose(fp_rom6);
	fclose(fp_rom7);
	fclose(fp_rom8);
	fclose(fp_rom9);
	fclose(fp_rom10);
	fclose(fp_rom11);
	fclose(fp_rom12);
	fclose(fp_rom13);
	fclose(fp_rom14);
	fclose(fp_rom15);

	root64K=_inverseRoot(FFT_SIZE);
	fp_rom1 = fopen("invROM1.rcf", "wt");
	fp_rom2 = fopen("invROM2.rcf", "wt");
	fp_rom3 = fopen("invROM3.rcf", "wt");
	fp_rom4 = fopen("invROM4.rcf", "wt");
	fp_rom5 = fopen("invROM5.rcf", "wt");
	fp_rom6 = fopen("invROM6.rcf", "wt");
	fp_rom7 = fopen("invROM7.rcf", "wt");
	fp_rom8 = fopen("invROM8.rcf", "wt");
	fp_rom9 = fopen("invROM9.rcf", "wt");
	fp_rom10 = fopen("invROM10.rcf", "wt");
	fp_rom11 = fopen("invROM11.rcf", "wt");
	fp_rom12 = fopen("invROM12.rcf", "wt");
	fp_rom13 = fopen("invROM13.rcf", "wt");
	fp_rom14 = fopen("invROM14.rcf", "wt");
	fp_rom15 = fopen("invROM15.rcf", "wt");
	for(int i=0; i<256; i++) {
		
		fprintf(fp_rom1,"%016llX\n", _power(root64K, i));
		fprintf(fp_rom2,"%016llX\n", _power(root64K, i*2));
		fprintf(fp_rom3,"%016llX\n", _power(root64K, i*3));
		fprintf(fp_rom4,"%016llX\n", _power(root64K, i*4));
		fprintf(fp_rom5,"%016llX\n", _power(root64K, i*5));
		fprintf(fp_rom6,"%016llX\n", _power(root64K, i*6));
		fprintf(fp_rom7,"%016llX\n", _power(root64K, i*7));
		fprintf(fp_rom8,"%016llX\n", _power(root64K, i*8));
		fprintf(fp_rom9,"%016llX\n", _power(root64K, i*9));
		fprintf(fp_rom10,"%016llX\n", _power(root64K, i*10));
		fprintf(fp_rom11,"%016llX\n", _power(root64K, i*11));
		fprintf(fp_rom12,"%016llX\n", _power(root64K, i*12));
		fprintf(fp_rom13,"%016llX\n", _power(root64K, i*13));
		fprintf(fp_rom14,"%016llX\n", _power(root64K, i*14));
		fprintf(fp_rom15,"%016llX\n", _power(root64K, i*15));
		/*
		format_binary(_power(root64K, i), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom1,"%d", s[j]);
		}
		fprintf(fp_rom1,"%d\n", s[63]);

		format_binary(_power(root64K, i*2), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom2,"%d", s[j]);
		}
		fprintf(fp_rom2,"%d\n", s[63]);

		format_binary(_power(root64K, i*3), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom3,"%d", s[j]);
		}
		fprintf(fp_rom3,"%d\n", s[63]);

		format_binary(_power(root64K, i*4), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom4,"%d", s[j]);
		}
		fprintf(fp_rom4,"%d\n", s[63]);

		format_binary(_power(root64K, i*5), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom5,"%d", s[j]);
		}
		fprintf(fp_rom5,"%d\n", s[63]);

		format_binary(_power(root64K, i*6), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom6,"%d", s[j]);
		}
		fprintf(fp_rom6,"%d\n", s[63]);

		format_binary(_power(root64K, i*7), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom7,"%d", s[j]);
		}
		fprintf(fp_rom7,"%d\n", s[63]);

		format_binary(_power(root64K, i*8), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom8,"%d", s[j]);
		}
		fprintf(fp_rom8,"%d\n", s[63]);

		format_binary(_power(root64K, i*9), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom9,"%d", s[j]);
		}
		fprintf(fp_rom9,"%d\n", s[63]);

		format_binary(_power(root64K, i*10), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom10,"%d", s[j]);
		}
		fprintf(fp_rom10,"%d\n", s[63]);

		format_binary(_power(root64K, i*11), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom11,"%d", s[j]);
		}
		fprintf(fp_rom11,"%d\n", s[63]);

		format_binary(_power(root64K, i*12), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom12,"%d", s[j]);
		}
		fprintf(fp_rom12,"%d\n", s[63]);

		format_binary(_power(root64K, i*13), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom13,"%d", s[j]);
		}
		fprintf(fp_rom13,"%d\n", s[63]);

		format_binary(_power(root64K, i*14), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom14,"%d", s[j]);
		}
		fprintf(fp_rom14,"%d\n", s[63]);

		format_binary(_power(root64K, i*15), s);
		for(int j=0; j<63; j++) {
			fprintf(fp_rom15,"%d", s[j]);
		}
		fprintf(fp_rom15,"%d\n", s[63]);
		*/

	}
	fclose(fp_rom1);
	fclose(fp_rom2);
	fclose(fp_rom3);
	fclose(fp_rom4);
	fclose(fp_rom5);
	fclose(fp_rom6);
	fclose(fp_rom7);
	fclose(fp_rom8);
	fclose(fp_rom9);
	fclose(fp_rom10);
	fclose(fp_rom11);
	fclose(fp_rom12);
	fclose(fp_rom13);
	fclose(fp_rom14);
	fclose(fp_rom15);

	largeFFT(a64, aFFT, FFT_SIZE, 0);
	largeFFT(b64, bFFT, FFT_SIZE, 0);
	//radix16FFT(a64, aFFT, FFT_SIZE, 4, 0);
	//radix16FFT(aFFT, bFFT, FFT_SIZE, 4, 1);
	//smallFFT16(a64, bFFT, FFT_SIZE, 0);
	//radix64FFT(a64, bFFT, FFT_SIZE, 1, 0);
	//radix64FFT(b64, bFFT, FFT_SIZE, 3, 0);

		for(int i=0; i<FFT_SIZE; i++) {
		d0 = i % 16;
		d1 = (i>>4) % 16;
		d2 = (i>>8) % 16;
		d3 = (i>>12) % 16;

		dataIndex = i >> 4;
		bankIndex = (d3+d2+d1+d0) % 16;
		switch(bankIndex){
			case 0:a64_s0[dataIndex] = aFFT[i];break; //a64[i]
			case 1:a64_s1[dataIndex] = aFFT[i];break; //a64[i]
			case 2:a64_s2[dataIndex] = aFFT[i];break; //a64[i]
			case 3:a64_s3[dataIndex] = aFFT[i];break; //a64[i]
			case 4:a64_s4[dataIndex] = aFFT[i];break; //a64[i]
			case 5:a64_s5[dataIndex] = aFFT[i];break; //a64[i]
			case 6:a64_s6[dataIndex] = aFFT[i];break; //a64[i]
			case 7:a64_s7[dataIndex] = aFFT[i];break; //a64[i]
			case 8:a64_s8[dataIndex] = aFFT[i];break; //a64[i]
			case 9:a64_s9[dataIndex] = aFFT[i];break;
			case 10:a64_s10[dataIndex] = aFFT[i];break;
			case 11:a64_s11[dataIndex] = aFFT[i];break;
			case 12:a64_s12[dataIndex] = aFFT[i];break;
			case 13:a64_s13[dataIndex] = aFFT[i];break;
			case 14:a64_s14[dataIndex] = aFFT[i];break;
			case 15:a64_s15[dataIndex] = aFFT[i];break;
		}
	}

	fp_data = fopen("InvdataIn.txt", "wt");
	for(int i=0; i<FFT_SIZE; i++) {
		fprintf(fp_data,"%016llX\n", aFFT[i]);
	}
	fclose(fp_data);

	fp_data0 = fopen("InvdataIn_0.txt", "wt");
	fp_data1 = fopen("InvdataIn_1.txt", "wt");
	fp_data2 = fopen("InvdataIn_2.txt", "wt");
	fp_data3 = fopen("InvdataIn_3.txt", "wt");
	fp_data4 = fopen("InvdataIn_4.txt", "wt");
	fp_data5 = fopen("InvdataIn_5.txt", "wt");
	fp_data6 = fopen("InvdataIn_6.txt", "wt");
	fp_data7 = fopen("InvdataIn_7.txt", "wt");
	fp_data8 = fopen("InvdataIn_8.txt", "wt");
	fp_data9 = fopen("InvdataIn_9.txt", "wt");
	fp_data10 = fopen("InvdataIn_10.txt", "wt");
	fp_data11 = fopen("InvdataIn_11.txt", "wt");
	fp_data12 = fopen("InvdataIn_12.txt", "wt");
	fp_data13 = fopen("InvdataIn_13.txt", "wt");
	fp_data14 = fopen("InvdataIn_14.txt", "wt");
	fp_data15 = fopen("InvdataIn_15.txt", "wt");
	for(int i=0; i<256; i++) {
		fprintf(fp_data0,"%016llX\n", a64_s0[i]);
		fprintf(fp_data1,"%016llX\n", a64_s1[i]);
		fprintf(fp_data2,"%016llX\n", a64_s2[i]);
		fprintf(fp_data3,"%016llX\n", a64_s3[i]);
		fprintf(fp_data4,"%016llX\n", a64_s4[i]);
		fprintf(fp_data5,"%016llX\n", a64_s5[i]);
		fprintf(fp_data6,"%016llX\n", a64_s6[i]);
		fprintf(fp_data7,"%016llX\n", a64_s7[i]);
		fprintf(fp_data8,"%016llX\n", a64_s8[i]);
		fprintf(fp_data9,"%016llX\n", a64_s9[i]);
		fprintf(fp_data10,"%016llX\n", a64_s10[i]);
		fprintf(fp_data11,"%016llX\n", a64_s11[i]);
		fprintf(fp_data12,"%016llX\n", a64_s12[i]);
		fprintf(fp_data13,"%016llX\n", a64_s13[i]);
		fprintf(fp_data14,"%016llX\n", a64_s14[i]);
		fprintf(fp_data15,"%016llX\n", a64_s15[i]);
	}
	fclose(fp_data0);
	fclose(fp_data1);
	fclose(fp_data2);
	fclose(fp_data3);
	fclose(fp_data4);
	fclose(fp_data5);
	fclose(fp_data6);
	fclose(fp_data7);
	fclose(fp_data8);
	fclose(fp_data9);
	fclose(fp_data10);
	fclose(fp_data11);
	fclose(fp_data12);
	fclose(fp_data13);
	fclose(fp_data14);
	fclose(fp_data15);

/*
	for(int i=0; i<FFT_SIZE; i++) {
		if(aFFT[i]!=bFFT[i]) {
			printf("test incorrect\n"); break;
		}
		//else
			//printf("test right");
	}
*/
	FILE *fp_dataOut = NULL;
	fp_dataOut = fopen("dataOut.txt", "wt");
	for(int i=0; i<FFT_SIZE; i++) {
		fprintf(fp_dataOut,"%016llX\n", aFFT[i]);
	}
	fclose(fp_dataOut);

	fp_dataOut = fopen("invFFT.txt", "wt");
	for(int i=0; i<FFT_SIZE; i++) {
		fprintf(fp_dataOut,"%016llX\n", bFFT[i]);
	} 
	fclose(fp_data);
	
	for(int i=0; i<FFT_SIZE; i++) {
		cFFT[i]=_multiply(aFFT[i], bFFT[i]);
	}

	largeFFT(cFFT, cFFT1, FFT_SIZE, 1);
	//radix64FFT(cFFT, cFFT1, FFT_SIZE, 3, 1);

	resolve_by3(cFFT1, FFT_SIZE);
	join_by3(cFFT1, FFT_SIZE, c, 2*words);
	
	fp_data = fopen("Mul_Result.txt", "wt");
	for(int i=0; i<2*words; i++) {
		fprintf(fp_data,"%X\n", c[i]);
	}
	fclose(fp_data);
		
	cZZ = ZZFromBytes((const unsigned char *)c, sizeof(uint32_t)*words*2);

	aZZ = ZZFromBytes((const unsigned char *)a, sizeof(uint32_t)*words);
	bZZ = ZZFromBytes((const unsigned char *)b, sizeof(uint32_t)*words);

	ZZ  cNtl;
	mul(cNtl, aZZ, bZZ);

	BytesFromZZ((unsigned char *)a, cZZ, sizeof(uint32_t)*words);
	BytesFromZZ((unsigned char *)b, cNtl, sizeof(uint32_t)*words);

	int i=0;

	if(compare(cZZ,cNtl)!=0) {
		printf("test is not correct!!!!!!\n");
		for(int i=0; i<words; i++){
			if(a[i]!=b[i]) {
				printf("a_fft[%d]=%x\n", i, a[i]);
				printf("a_ntl[%d]=%x\n", i, b[i]);
			}
		}
		//printf("a_fft[%d]=%x\n", i, a[i]);
		//printf("a_ntl[%d]=%x\n", i, b[i]);
		//printf("a_fft[%d]=%x\n", i+1, a[i+1]);
		//printf("a_ntl[%d]=%x\n", i+1, b[i+1]);
	}
	else printf("test is correct.\n");
	


	free(a);
	free(b);
	free(c);
	free(aFFT);
	free(bFFT);
	free(cFFT);
	free(a64);
	free(b64);
		
	getchar();
}
