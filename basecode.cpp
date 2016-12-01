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
#include <assert.h>
#include <vector>
using namespace std;
using namespace NTL;

#define _N        0xFFFFFFFF00000001ull
#define _ROOT     0xE9653C8DEFA860A9ull    /* 3*2^32th root of unity */
#define _CRITICAL 0xFFFFFFFFul
#define _N2       0x0000000000010001ull //65537
#define _ROOT2 	  0x3ull

void output(uint64_t* x, int length){
	for(int i = 0; i < length; i++){
		cout << x[i] << ", ";
	}
	cout << endl << endl;
	return;
}

void output(uint32_t* x, int length){
	for(int i = 0; i < length; i++){
		cout << x[i] << ", ";
	}
	cout << endl << endl;
	return;
}
uint64_t _add(uint64_t a, uint64_t b) {
	uint64_t sum = a + b;
  	if(sum < a){
  	      	if(sum >= _N)
  	         	sum = sum - _N - _N;
  	       	else
  	         	sum = sum - _N;
  	      }
  	else if(sum >= a && sum >= _N)
  	      	sum = sum - _N;

  	return sum;
}

uint64_t _add2(uint64_t a, uint64_t b) {
//	uint64_t sum = a + b;
//  if(sum < a){
//		if(sum >= _N2)
//	   	sum = sum - _N2 - _N2;
//	 	else
//	   	sum = sum - _N2;
//	}
//  else if(sum >= a && sum >= _N2)
//		sum = sum - _N2;
//
//  return sum;
//
	assert(a < _N2 * _N2 && b < _N2 * _N2);
	uint64_t sum = a + b;
	while(sum >= _N2){
		sum -= _N2;
	}
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
uint64_t _subtract2(uint64_t a, uint64_t b) {
  uint64_t sum;

  if(b>=_N2)
    b-=_N2;
  if(a>=b)
    sum=a-b;
  else
    sum=a-b+_N2;
  
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

uint64_t _multiply2(uint64_t a, uint64_t b) {
	assert(a < (_N2*_N2)  && b <  (_N2*_N2));
	uint64_t ab = a * b;
	uint64_t ab2 = (ab & 0xffff00000000ull) >> 32;
	uint64_t ab1 = (ab & 0x0000ffff0000ull) >> 16;
	uint64_t ab0 = (ab & 0x00000000ffffull);

	uint64_t cSubB;
	uint64_t res;
	cSubB = _subtract2(ab0, ab1);
	res = _add2(cSubB, ab2);	
	return res;

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
uint64_t _power2(uint64_t a, uint64_t k) {
  uint64_t current=1, square=a;

  while(k>0) {
    if((k&1) != 0)
      current = _multiply2(current, square);
    square = _multiply2(square, square);
    k = k >> 1;
  }
  
  return current;
}

uint64_t _normalize(uint64_t a) {
  if(a>=_N)
    return a-_N;
  return a;
}

uint64_t _normalize2(uint64_t a) {
  while(a>=_N2){
		a  = a - _N2;	
	}
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
uint64_t _inverse2(uint64_t a) {
  int64_t s, t;

  extendedGCD(_N2, a, &s, &t);
  if(t<0)
    return t+_N2;
  return t;
}

uint64_t _root(uint64_t size) {
  uint64_t k=(3ull<<32)/size;
  return _normalize(_power(_ROOT, k));
}

uint64_t _root2(uint64_t size) {
  uint64_t k=(65536)/size;
  return _normalize2(_power2(_ROOT2, k));
}

uint64_t _inverseRoot(uint64_t size) {
  uint64_t k=(3ull<<32)/size;
  return _normalize(_power(_ROOT, (3ull<<32)-k));
}

uint64_t _inverseRoot2(uint64_t size) {
  uint64_t k=(65536)/size;
  return _normalize2(_power2(_ROOT2, (65536)-k));
}

void transpose(uint64_t *x, uint64_t *X, uint32_t size, uint32_t xLength) {
  int group, index;

  for(group=0;group<size/xLength;group++)
    for(index=0;index<xLength;index++)
      X[group*xLength+index]=x[index*size/xLength+group];
}
void transpose64K1(int size) {
  for(int index = 0; index < size; index++){
  	cout << index << " -> " << (index % 4096) * 16 + index/4096 << endl;
  }

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
	//	X[i] = _normalize(value); 
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
void smallFFT_debug(uint64_t *x, uint64_t* temp, uint64_t *X, uint32_t size, int inverse) {
	int i, j;
	uint64_t total, r, inv=1, value;
	if(inverse) {
		r=_inverseRoot(size);
		inv=_inverse(size);
	}
	else
		r=_root(size);
	cout << "root64 is: " << r << endl;
	//printf("root8=%d\n", r);
	for(i=0;i<size;i++) {
		total=0;
		for(j=0;j<size;j++) {
			value=_multiply(_power(r, i*j), x[j]);
			total=_add(total, value);
			temp[i * 64 + j] = value; 
		}
		X[i]=_normalize(_multiply(total, inv));
	}
	//printf("smallFFT:X[0]=%x\n",X[0]);
}
void smallFFT2(uint64_t *x, uint64_t *X, uint32_t size, int inverse) {
  int i, j;
  uint64_t total, r, inv=1, value;
  if(inverse) {
    r=_inverseRoot2(size);
    inv=_inverse2(size);
  }
  else
    r=_root2(size);
    //printf("root8=%d\n", r);
    for(i=0;i<size;i++) {
    total=0;
    for(j=0;j<size;j++) {
      value=_multiply2(_power2(r, i*j), x[j]);
      total=_add2(total, value);
    }
    X[i]=_normalize2(_multiply2(total, inv));
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
	int count = 1, fft, i, j, k, l, m;
	uint64_t r, inv=1, omega, currentOmega, c0, c1;
	uint64_t *buffer, *to, *from, *swap, *smallX, *smallY;
	buffer = (uint64_t *)malloc(sizeof(uint64_t)*size);

	if(size%3==0) {
		size = size / 3;
		count = 3;
	}
	if(inverse) {
		r = _inverseRoot(size);
		inv = _inverse(size);
	}
	else
		r = _root(size); //
	

	for(fft = 0;fft < count; fft++) {
		for(i = 0;i < size; i++)
			buffer[i] = x[count * i + fft];
		
		from = buffer;
		to = X + fft * size;
		l = size / 2;
		m = 1;
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

void largeFFT2(uint64_t *x, uint64_t *X, uint32_t size, int inverse) {
  int      count=1, fft, i, j, k, l, m;
  uint64_t r, inv=1, omega, currentOmega, c0, c1;
  uint64_t *buffer, *to, *from, *swap, *smallX, *smallY;

  if(size%3==0) {
    size=size/3;
    count=3;
  }
  if(inverse) {
    r=_inverseRoot2(size);
    inv=_inverse2(size);
  }
  else
    r=_root2(size); //
		
  buffer=(uint64_t *)malloc(sizeof(uint64_t)*size);

  for(fft=0;fft<count;fft++) {
    for(i=0;i<size;i++)
      buffer[i]=x[count*i+fft];
    from=buffer;
		to=X+fft*size;
		l=size/2;
		m=1;
		while(l>=1) {
  		omega=_power2(r, size/(2*l));
  		currentOmega=1;
			for(j=0;j<l;j++) {
				for(k=0;k<m;k++) {
					c0=from[k+j*m];
					c1=from[k+j*m+l*m];
					to[k+2*j*m]=_add2(c0, c1);
					to[k+2*j*m+m]=_multiply2(currentOmega, _subtract2(c0, c1));
				}
				currentOmega=_multiply2(currentOmega, omega);
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
		r=_root2(size*count);
		smallX=(uint64_t *)malloc(sizeof(uint64_t)*count);
		smallY=(uint64_t *)malloc(sizeof(uint64_t)*count);
		for(fft=0;fft<size;fft++) {
			smallX[0]=X[fft];
			smallX[1]=_multiply2(_power2(r, fft), X[fft+size]);
			smallX[2]=_multiply(_power2(r, fft+fft), X[fft+size+size]);
			smallFFT2(smallX, smallY, 3, inverse);
			X[fft]=smallY[0];
			X[fft+size]=smallY[1];
			X[fft+size+size]=smallY[2];
		}
		free(smallX);
		free(smallY);
	}

	for(i=0;i<size*count;i++)
		X[i]=_normalize2(_multiply2(X[i], inv));
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

		// cout << endl << endl << "for the " << i << "th level" << endl;
		for(j=0;j<l;j++) {
			if(inverse)
				r = _inverseRoot(64*l);
			else
				r = _root(64*l);
			
				for(k=0;k<m;k++) {
					// cout << endl << endl;
					for(int t=0; t<64; t++){
						//if(t == 0 || t == 1 || t == 2|| t == 63 || t == 62)
						// cout << "	x[" << k + j*m + t*l*m << "] -> x[" << t << "]" << endl; 
          				smallX[t]=from[k + j*m + t*l*m];
					}
					// cout << endl << "	64 FFT " << endl << endl;
					smallFFT(smallX, smallY, 64, inverse);

					for(int t=0; t<64; t++) {
						if(inverse)
							to[k+64*j*m+t*m]=_multiply(smallY[t], _power(r, t * j));
						else{
							//if(t == 0 || t == 1 || t == 2|| t == 63 || t == 62)
								// cout << "	X[" << t << "]" << "*W" << 64*l << "^" << t * j << " -> " << k+64*j*m+t*m << endl;
							to[k+64*j*m+t*m]=_multiply(smallY[t], _power(r, t * j));
						}
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

	free(buffer);
}


void largeFFT64(uint64_t *x, uint64_t *X, uint64_t size, int inverse) {
	uint64_t l = size / 4096;
	uint64_t* bufferIn4K =(uint64_t *)malloc(sizeof(uint64_t)*4096);
	uint64_t* bufferOut4K =(uint64_t *)malloc(sizeof(uint64_t)*4096);	
	uint64_t* buffer = (uint64_t *)malloc(sizeof(uint64_t)*size);	
	uint64_t smallX[16], smallY[16];

	uint64_t w;

	for(int i = 0; i < l; i++){
		// cout << endl << endl<< "for the " << i + 1 << "th level" << endl;
		if(inverse){
			w = _inverseRoot(size);
		}
		else{
			w = _root(size);
		}

		for(int j = 0; j < 4096; j++){
			// if(j == 0 || j == 1 || j == 2 || j == 4095 || j == 4094){
			// 	cout << "	x[" << j * 16 + i << "] -> " << "x[" << j << "]" << endl;
			// }
			bufferIn4K[j] = x[j * 16 + i];
		}
		// cout << endl << "	4K FFT" << endl << endl;
		//largeFFT(bufferIn4K, bufferOut4K, 4096, inverse);
		radix64FFT(bufferIn4K, bufferOut4K, 4096, 2, 0);

		for(int j = 0; j < 4096; j++){
			// if(j == 0 || j == 1 || j == 2 || j == 4095 || j == 4094){
			// 	cout << "	X[" << j << "] * W" << size << "^" << i * j << " -> buffer[" << i * 4096 + j<<"]" << endl;
			// }
			buffer[i * 4096 + j] = _multiply(bufferOut4K[j], _power(w, j * i));
		}
	}
	cout << "**********************************" << endl;
	output(buffer, 64);
	cout << "**********************************" << endl;

	for(int i = 0; i < size/16; i++){
		// cout << endl << endl<< "for the " << i + 1 << "th 16" << endl;

		for(int j = 0; j < 16; j++){
			cout << "	buffer[" << j * 4096 + i << "] -> x[" << j << "]" << endl;
			smallX[j] = buffer[j * 4096 + i];
		}

		// cout << endl << "	16 FFT" << endl << endl;
		smallFFT(smallX, smallY, 16, inverse);

		for(int j = 0; j < 16; j++){
			// cout << "	X[" << j << "] -> result[" << j * 4096 + i << "]" << endl;
			X[j * 4096 + i] = smallY[j];
		}
	}


	free(bufferIn4K); free(bufferOut4K); free(buffer);
	bufferIn4K = NULL; bufferOut4K = NULL; buffer = NULL;
	
	return;
}
void smallFFT16(uint64_t *x, uint64_t *X, uint32_t size, int inverse) {
	int i, j;
	uint64_t total, r, inv=1, value;

	if(inverse) {
		r=_inverseRoot(size);
		inv=_inverse(size);
	}
	else
		r=_root(size);

	for(i=0;i<size;i++) {
		total=0;
		for(j=0;j<size;j++) {
			value=_multiply(_power(r, i*j), x[j]);
			total=_add(total, value);
		}
		X[i]=_normalize(_multiply(total, inv));
	}
}

void radix416(uint64_t *x, uint64_t *X, uint32_t size, int inverse){
	int l, j, k;
	uint64_t *buffer, *to, *from, *swap;
	uint64_t smallX[4], smallY[4], r;
	
	buffer=(uint64_t *)malloc(sizeof(uint64_t)*size);
	for(int i = 0; i < size; i++){
		buffer[i]=x[i];
	}
  	

	from = buffer;
	to = X;
	l = size/4;
	int m = 1;

	for(int i = 0; i < 2; i++){
		cout << "for the " << i+1 << " th level: " << endl;
		for(j = 0; j < l; j++){
			if(inverse)
	    		r = _power(_inverseRoot(4 * l), j);
			else
				r = _power(_root(4 * l), j);

			for(k = 0; k < m; k++) {
				cout /*<< "	in " << k << " of m:"*/ << endl;
				cout << " transposing as follow: " << endl;
				for(int t = 0; t < 4; t++){
					smallX[t] = from[k + j*m + t*l*m];
						cout << "		x[" << k + j*m + t*l*m << "] -> x[" << t << "]" << endl;
				}

				smallFFT(smallX, smallY, 4, inverse);
					cout << "......" << endl << "		x[i]" << " 16 point DFT X[i]" << endl << "......" << endl;
				for(int t=0; t<4; t++) {
					if(inverse)
						to[k + 4*j*m + t*m] = _multiply(smallY[t], _power(r, t));
					else{
						if(i == 0)
							to[k + 4*j*m + t*m] = _multiply(smallY[t], _power(r, t)); //
						else
							to[k + 4*j*m + t*m] = smallY[t];
						
						if (i == 0)
							cout << "		" << "X[" << t << "]" << " * W" << 4*l << "^" << t << " -> " << k + 4*j*m + t*m << endl;
						else
							cout << "		" << "X[" << t << "]" << " -> " << k + 4*j*m + t*m << endl;
					}
				}
			} 
		}
		swap = from;
		from = to;
		to = swap;
		l = l/4;
		m = m*4;

	}
	for(int i=0;i<size;i++)
		X[i] = from[i];
	
	free(buffer);
}

void radix16FFT(uint64_t *x, uint64_t *X, uint32_t size, uint32_t level, int inverse) {
   
  uint32_t l, m, j, k, base, current, offset;
  uint64_t inv = 1, omega, currentOmega, c0, c1;
  uint64_t *buffer, *to, *from, *swap;
  uint64_t r;
  uint64_t smallX[16], smallY[16];
  uint64_t temp0 = 0;
     
  buffer=(uint64_t *)malloc(sizeof(uint64_t)*size);

	for(int t = 0; t < size; t++)
		buffer[t]=x[t];
	from=buffer;
	to=X;
	l = size/16;
    m=1;
   
	for(int i = 0; i < level; i++) {
		cout << "for the " << i << " th level: " << endl;
		for(j = 0; j < l; j++) {
				cout /*<< "	in " << j << "th 16:"*/ << endl;  
			if(inverse)
        		r = _power(_inverseRoot(16 * l), j);
			else
				r = _power(_root(16 * l), j);
			
			for(k = 0; k < m; k++) {
				cout /*<< "	in " << k << " of m:"*/ << endl;
				cout << " transposing as follow: " << endl;
				for(int t = 0; t < 16; t++){
					smallX[t] = from[k + j*m + t*l*m];
						cout << "		x[" << k + j*m + t*l*m << "] -> x[" << t << "]" << endl;
				}

				smallFFT16(smallX, smallY, 16, inverse);
					cout << "......" << endl << "		x[i]" << " 16 point DFT X[i]" << endl << "......" << endl;
				for(int t=0; t<16; t++) {
					if(inverse)
						to[k + 16*j*m + t*m] = _multiply(smallY[t], _power(r, t));
					else{
						if(i == 0)
							to[k + 16*j*m + t*m] = _multiply(smallY[t], _power(r, t)); //
						else{
							to[k + 16*j*m + t*m] = smallY[t];
						}
							cout << "		" << "X[" << t << "]" << " * W" << 16*l << "^" << t * j << " -> " << k + 16*j*m + t*m << endl;

					}
				}
			} 
		}
		swap = from;
		from = to;
		to = swap;
		l = l/16;
		m = m*16;
	}
	for(int i=0;i<size;i++)
		X[i] = from[i];
	
	free(buffer);
}

int format_binary(uint64_t x, uint32_t *s)
{
	for (int z = 0; z < 64; z++) {
		s[63-z] = ((x>>z) & 0x0000000000000001) ? 1 : 0;
	} 
	return 0;
}

ZZ convToZZ(int size, uint64_t* data, int significantBytes){
  ZZ res;
	char* buff = (char*)malloc(sizeof(char) * size * significantBytes);

	for(int i = 0; i < size; i++){
		memcpy(buff + i*significantBytes, data + i, sizeof(char) * significantBytes);
	}   
	res = ZZFromBytes((const unsigned char*)buff, size * significantBytes);

	free(buff); 
	buff = NULL;
	return res;
}

vector<ZZ> CRTmul(uint64_t* x, uint64_t*y, int size, int fft_size){
	assert(size * 2 <= fft_size);
	uint64_t* X1 = (uint64_t*)malloc(sizeof(uint64_t) * fft_size);
	uint64_t* X2 = (uint64_t*)malloc(sizeof(uint64_t) * fft_size);
	uint64_t* Y1 = (uint64_t*)malloc(sizeof(uint64_t) * fft_size);
	uint64_t* Y2 = (uint64_t*)malloc(sizeof(uint64_t) * fft_size);
	uint64_t* C1 = (uint64_t*)malloc(sizeof(uint64_t) * fft_size);
	uint64_t* C2 = (uint64_t*)malloc(sizeof(uint64_t) * fft_size);
	uint64_t* c1 = (uint64_t*)malloc(sizeof(uint64_t) * fft_size);
	uint64_t* c2 = (uint64_t*)malloc(sizeof(uint64_t) * fft_size);
	
	largeFFT(x, X1, fft_size, 0);
	largeFFT(y, Y1, fft_size, 0);
	largeFFT2(x, X2, fft_size, 0);
	largeFFT2(y, Y2, fft_size, 0);

	for(int i = 0; i < fft_size; i++){
		C1[i] = _multiply(X1[i], Y1[i]);
		C2[i] = _multiply2(X2[i], Y2[i]);
	}		
	
	largeFFT(C1, c1, fft_size, 1);
	largeFFT2(C2, c2, fft_size, 1);

	uint64_t carry = 0;
	for(int i = 0; i < fft_size; i++){
		c1[i] += carry;
		uint64_t current = c1[i];
		carry = current >> 32;
		c1[i] = c1[i] & 0xFFFFFFFF;
	}
	carry = 0;
	for(int i = 0; i < fft_size; i++){
		c2[i] += carry;
		uint64_t current = c2[i];
		carry = current >> 32;
		c2[i] = c2[i] & 0xFFFFFFFF;
	}
	for(int i =0; i < 10; i++){
		cout << c1[i] << ", " << c2[i] << endl;
	
	}

	ZZ N1; uint64_t N1t = _N; 
	N1 = convToZZ(1, &N1t, 8);
	ZZ N2;N2 = _N2;
	ZZ N1N2;
	ZZ term1, term2;
	mul(N1N2, N1, N2);

	vector<ZZ> result(fft_size);


	//for(int i = 0; i < fft_size; i++){
	//	uint64_t tSubS =  _subtract2(c2[i], c1[i]);
	// 	uint64_t inverse = _inverse2(_N);
	// 	term1  = inverse * tSubS;
	// 	term1 = term1 * N1;
	// 	term2 = c1[i]; 
	// 	result[i]  = (term1 + term2) % N1N2;
	//}
	
	uint64_t y1 = _inverse(_N2);
	cout << "the inverse of _N2 is: " << y1 << endl;
	uint64_t y2 = _inverse2(_N);
	cout << "the inverse of _N1 is: " << y2 << endl;
	ZZ y1ZZ, y2ZZ;
	y1ZZ = y1; y2ZZ = y2;
	for(int i = 0; i < fft_size; i++){
		term1 = convToZZ(1, &c1[i], 4);
		term2 = convToZZ(1, &c2[i], 4);
		term1 *= y1ZZ;
		term2 *= y2ZZ;
		term1 *= N2;
		term2 *= N1;
		result[i]  = (term1 + term2) % N1N2;
		if(result[i] < 0){
			cout << "here!!!!!!!!!!!!" << endl;
			cout << N1N2 << endl;
			result[i] += N1N2;
		}
	}

	free(X1);
	free(X2);
	free(Y1);
	free(Y2);
	free(C1);
	free(C2);
	free(c1);
	free(c2);

	return result;
}

bool barrettPreComputation(uint64_t& n, uint64_t& u, uint64_t m){
  n = NextPowerOfTwo(m);
  
  if(n > 31){
    return false;
  }

  long powerOfTwo = 1ull << (2*n);
  u = powerOfTwo / m;

  return true;
}

void barrettMul(uint32_t& a, uint32_t& b, uint64_t& result, uint64_t n, uint64_t u, uint64_t m){
  uint64_t ab = a * b;
  uint64_t r = ab - (ab >> n) * (u >> n) * m;
  uint64_t count = 0;
  while(r >= m){
    r -= m;
    count++;
  }
  cout << count << endl;
  result = r;
  return;
}


void twiddleGen(int size,  uint64_t* input){
	uint64_t r = _root(size);
	int len = sqrt(size);
	cout << len << endl;
	for(int i = 0; i < len; i++){
		for(int j = 0; j < len; j++){
			input[i * 64 + j] = _power(r, i * j);	
		}
	}
	return;
}

void twiddleGen64K(int size,  uint64_t* input){
	uint64_t r = _root(size);
	int len = size / 4096;
	
	for(int i = 0; i < len; i++){
		for(int j = 0; j < 4096; j++){
			input[i * 4096 + j] = i * j;	
		}
	}
	return;
}

int main(){

	uint64_t fftSize = 4096;

	uint64_t* x = (uint64_t*)malloc(sizeof(uint64_t) * fftSize);
	uint64_t* X = (uint64_t*)malloc(sizeof(uint64_t) * fftSize);
	
	for(int i = 0; i < fftSize; i++){
		x[i] = i;
	}

	largeFFT(x,X,fftSize,0);
	output(X,fftSize);
	free(x);free(X);
	x = NULL; X = NULL; 
	return 0;
}
