#include <cstring>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <bitset>
#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <sys/time.h>


using namespace std;
using namespace NTL;

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

  upper = (albl>>32) + albh + ahbl;
  if(upper<ahbl)
    carry=0x100000000ull;
  upper = (upper>>32) + ahbh + carry;

  uu = upper>>32;
  ul = upper;

  if(ul==0)
    upper = _N-uu;
  else
    upper = (upper<<32)-uu-ul;

  return _add(a * b, upper);
}

uint64_t _power(uint64_t a, uint64_t k) {
  uint64_t current=1, square=a;

  while(k>0) {
    if((k&1) != 0)
      current = _multiply(current, square);
    square = _multiply(square, square);
    k = k >> 1;
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
    int i, j;
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
  int i, j;
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
  int index, byteCount = 0, current = 0;
  uint32_t value = 0;
  for(index = 0;index < XLength; index++){ 
    if(byteCount == 0 && current < xLength) {
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
  uint64_t carry=0, current;
  int      index;
  for(index=0; index<length; index++) {
    current=x[index];
    carry+=(current & 0xFFFFFF);
    x[index]=carry & 0xFFFFFF;
    carry=(carry>>24) + (current>>24);
  }
  if(carry>0) {
    printf("resolve_by3 failed! abort!\n");
    exit(1);
  }
}

void join_by2(uint64_t *x, uint32_t xLength, uint32_t *X, uint32_t XLength) {
  int      byteCount=0, current=0, index;
  uint64_t value=0;
  value=0;
  byteCount=0;
  current=0;
  for(index=0;index<XLength;index++) {
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

  //printf("m=%d\n",m);

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
  //printf("smallFFT:inv=%016llX\n",inv);
}

// void radix16FFT(uint64_t *x, uint64_t *X, uint32_t size, uint32_t level, int inverse) {
   
//   uint32_t l, m, j, k, base, current, offset;
//   uint64_t inv=1, omega, currentOmega, c0, c1;
//   uint64_t *buffer, *to, *from, *swap;
//   uint64_t r;
//   uint64_t smallX[16], smallY[16];
//   uint64_t temp0=0;
     
//   buffer=(uint64_t *)malloc(sizeof(uint64_t)*size);

//   for(int t=0;t<size;t++)
//     buffer[t]=x[t];
//     from=buffer;
//     to=X;
//     l=size/16;
//     m=1;
   
//   for(int i=0; i<level; i++) {
//     FILE *fp_dataOut = NULL;
//     fp_dataOut = fopen("dataStageOut.txt", "wt");
//     for(j=0;j<l;j++) {
//       if(inverse)
//         r=_power(_inverseRoot(16*l), j);
//       else
//         r=_power(_root(16*l), j);
// 			for(k=0;k<m;k++) {
// 				for(int t=0; t<16; t++)
// 				  smallX[t]=from[k+j*m+t*l*m];
//   				//for(int t=0; t<16; t++) {
//   				//	printf("index=%d, X = %016llX\n", k+j*m+t*l*m, smallX[t]);
//   				//}
// 				smallFFT16(smallX, smallY, 16, inverse);
// 				for(int t=0; t<16; t++) {
// 					if(inverse)
// 						to[k+16*j*m+t*m]=_multiply(smallY[t], _power(r, t));
// 					else
// 						to[k+16*j*m+t*m]=_multiply(smallY[t], _power(r, t)); //
// 				}
// 				for(int t=0; t<16; t++) {
// 					//printf("index=%d, smallY=%016llX Pow=%016llX Y = %016llX\n", k+j*m+t*l*m, smallY[t], _power(r, t), _multiply(smallY[t], _power(r, t)));
// 					//fprintf(fp_dataOut, "index=%d, smallX=%016llX smallY=%016llX Pow=%016llX Y = %016llX\n", k+j*m+t*l*m,smallX[t], smallY[t], _power(r, t), _multiply(smallY[t], _power(r, t)));
// 				}
// 			} 
// 		}


//     //for(int i=0; i<size; i++) {
//     //    fprintf(fp_dataOut,"%016llX\n", to[i]);
//     //}
//     fclose(fp_dataOut);
// 		swap=from;
// 		from=to;
// 		to=swap;
// 		l=l/16;
//    m=m*16;
//   }
//   for(int i=0;i<size;i++)
//     X[i]=from[i];
// }

int format_binary(uint64_t x, uint32_t *s)
{
  for (int z = 0; z < 64; z++) {
    s[63-z] = ((x>>z) & 0x0000000000000001) ? 1 : 0;
  } 
  return 0;
}

ZZ convToZZ(int size, uint64_t* data){
    ZZ res;
    char* buff = (char*)malloc(sizeof(char) * size * 3);

    for(int i = 0; i < size; i++){
      memcpy(buff + i*3, data + i, sizeof(char) * 3);
    }   
    res = ZZFromBytes((const unsigned char*)buff, size * 3);

    free(buff); 
    buff = NULL;
    return res;
}


int main(){
	int fftSize=16;

	uint64_t* x=(uint64_t*)malloc(sizeof(uint64_t)*fftSize);
	uint64_t* y=(uint64_t*)malloc(sizeof(uint64_t)*fftSize);
	
	for(int i=0;i<4;i++){
				x[i]=11;
	}
	for(int i=4;i<fftSize;i++){
				x[i]=33;
	}
  largeFFT(x,y,fftSize,0);
	cout<<"the radixj16 output are:"<<endl;
	for(int i=0;i<fftSize;i++){
	cout<<y[i]<<endl;}
	cout<<endl<<endl;
	free(x);free(y);
	x=NULL;y=NULL;
  return 0;
}
