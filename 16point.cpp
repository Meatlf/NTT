/*功能：公式（1）的实现*/
#include <NTL/ZZ.h>
using namespace NTL;
#include <math.h>
#include<iostream>
using namespace std;


int main()

{	/*定义长度为16的矢量x和y*/
	Vec<ZZ> x;
	Vec<ZZ> y;
	x.SetLength(16);
	y.SetLength(16);
	
	/*输入矢量x*/
	cin>>x ;
	
	/*定义计算p需要的变量*/
	ZZ  a, b, p;
	ZZ t;
	ZZ w;

	/*计算p*/
	power2(a, 64);
	power2(b, 32);
	p = a - b + 1;
	
	for (long k = 0; k < 16; k++)
	{
		y[k] = 0;
		for (long n = 0; n < 16; n++)
		{
			
			power2(t, (12 * n*k) % 192); w = x[n] * t;
			y[k] = y[k] +w ;
		}
		    y[k] = y[k] % p;
	}
	cout << y << '\n';
	return 0;
}
