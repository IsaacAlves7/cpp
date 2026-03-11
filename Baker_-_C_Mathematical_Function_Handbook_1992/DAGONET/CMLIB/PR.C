/*  powers and roots
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

power		x^n for x real general n (integers special)
root		x^(1/n)  for integral n
square_rt	square root
cube_rt		cube root
Euclidd		Euclidean distance. hypot() on some systems. 
		use in a "bulletproofed" form of cabs()
*/

#define tol 1.e-8
#include "cmlib.h"
#include "protom.h"

double power( double x, double n)
	{
	int k,i;double pow,p;
	if(n==0.)return 1.0;/* return x^0=1 (even if x=0)*/
	if(x==1.)return 1.;
	if(x==0.)
		{if(n>0.)return 0.;
		return errorcode;
		}
	k=n;
	/*non-integral*/
	if( ((double)k-n)!=0.)
		{
		if(x<=0.)return errorcode;
		return expon(n*ln(x));
		}
	pow=1.;
	p=x;if(x<0)p=-p;
	if(k<0){k=-k;i=1;}
	else i=0;
	while(k)
		{
		if( 1 & k)pow*=p;
		p*=p;
		k>>=1;		
		}
	if(x<0 && (((int)n)%2))
		{
		if(i) return -1./pow;
		return -pow;
		}
	if(i) return 1./pow;
	return pow;
	}


double root(xx,ni) double xx;int ni;
	{double mult,a,z,x,b,c,zo;int d,n;
	x=xx;mult=1.;n=abs(ni);
	if(n==0)return 1.0;/* return x^0=1 (even if x=0)*/
	if(x==0.)return 0.0;
	if(x==1.)return 1.;
	if(x<0 && !(n%2))return errorcode;
	if(x<0){x=-x;mult=-1.;}
	a=1./n;b=1.-a;c=a*x;d=n-1;
	z=zo=1.;
	for(;;)
		{
		z= b*z+c/power(z,d);
		if(abs(z-zo)<tol)
			{if(ni>0)return z*mult;
			return mult/z;
			}
		zo=z;
		}
	}


double square_rt(double x)
	{double /*y,*/z,old;
	if(x<0.)return errorcode;
	if(x==0.)return 0.;
	if(x==1.)return 1.;
	z=1.;/*initial guess*/
	infinite_loop
		{
		old=z;
		z= .5*(z+x/z);
		if( abs(old-z)<tol)break;
		}
	return z;
	}

double cube_rt(double x)
	{double /*y,*/z,old,c;
	if(x==0.)return 0.;
	if(x==1.)return 1.;
	c=1./3.;
	z=(x>0.)?1.:-1.;/*initial guess*/
	infinite_loop
		{
		old=z;
		z= (2.*z+x/(z*z))*c;
		if( abs(old-z)<tol)break;
		}
	return z;
	}

double Euclidd(xin,yin)double xin,yin;
	{
	double x,y,hold;
	x=abs(xin);y=abs(yin);	
	if(x<y){hold=x;x=y;y=hold;}
	/* x largest mag. of xin or yin */
	if(x==0.)return 0.;
	if(y==0.) return x;
	hold= y/x;
	return x*sqrt(1.+hold*hold);
	}
