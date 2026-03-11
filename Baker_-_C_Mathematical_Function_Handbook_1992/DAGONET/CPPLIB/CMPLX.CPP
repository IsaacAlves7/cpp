/*
Complex C++ Utilities

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

*/
#include "Cmatrix.hpp"

#define maxterm 100
#define abstol 1.e-10
#define reltol 1.e-5
#define zerotol 1.e-14
#define errorcode -1.e60
#define pi 3.141592653589793238462643383279
#define lim 15.
#define max(a,b) ((a)>(b)?(a):(b))
#define infinite_loop for(;;)
#define rabs(x) ((x)<0.? -(x):(x))
#define ABS(X) rabs(X)
#define topexp 350.

void complex::print( char *ahead,char *behind)
	{char *between="";
	if(y>=0.)between="+";
	printf("%s %e%s%e i %s",ahead,x,between,y,behind);}

double complex::abs()
	{ return sqrt(x*x+y*y);}

complex complex::cexp()
	{double scale;
	scale= x;
	if( x> topexp) return complex( errorcode,0.);
	if( x<-topexp) return complex( 0.,0.);
	scale= exp((scale));
	return complex(scale*cos(y),scale*sin(y));
	}
complex complex::clog()
	{double mant,arg,mag;mag=(*this).abs();
	if( mag < zerotol ) return complex(errorcode,0.);
	mant = log(mag);
	arg= atan2(y,x);
	return complex(mant,arg);
	}

complex complex::operator^(double expon)
	{
	complex z;
	z= (*this).clog( ) * expon;
	return z.cexp();
	}

complex complex::operator^(complex expon)
	{
	complex z;
	z= (*this).clog( ) * expon;
	return z.cexp();
	}


complex complex::csin()
{double z,zi,sinh,cosh,real;real=x;
if( y> topexp) return complex( errorcode,0.);
if( y<-topexp) return complex( 0.,0.);
z=exp( y); zi=1./z; cosh=.5*(z+zi);sinh=.5*(z-zi);
//printf(" csin real imag arg %e %e\n",real,z);
return complex(cosh*sin(real),cos(real)*sinh);
}

complex complex::ccos()
{double z,zi,sinh,cosh,real;real=x;
if( y> topexp) return complex( errorcode,0.);
if( y<-topexp) return complex( 0.,0.);
z=exp(y); zi=1./z; cosh=.5*(z+zi);sinh=.5*(z-zi);
return complex(cosh*cos(real),-sin(real)*sinh);
}

complex tangent ( complex& x)
{
return x.csin()/x.ccos();
}

complex tanh(complex& x)
{complex i(0.,1.);
return -i*tangent( x*i);
}
complex sinh(complex& x)
{complex i(0.,1.);
return  -i* ( x*i).csin();
}

complex asin(complex&x)
{complex i(0.,1.);
return      -i * (i*x+ ((1.-x*x)^.5)).clog();
}
complex acos(complex&x)
{complex i(0.,1.);
return      -i * (x+ i*((1.-x*x)^.5)).clog();
}

complex arctan( complex& x)
{complex i(0.,1.),a;
 a= ( i+x )/(i-x);
 return .5 * i * a.clog();
}


complex cgamma( complex& xin)
{complex xi,factor(1.,0.),g,lg,x;
//reflect
x=xin;
if(x.real()<0.)
	{return pi/(cgamma(1.-x)*(pi*x).csin());
	}
while(x.abs()< lim){factor=factor/x;x++;}xi=1./x;
g= factor*((-x).cexp())*((x)^(x-.5))*sqrt(2.*pi)*
(1.+xi*(1./(12.)+xi*(1./(288)-xi*(139./(51840)+571./(2488320.*x)))));
lg=g.clog();
//lg.print(" log gamma=","\n");
return g;
}


double loggam(x) double x;
{
int i;
double z,tmp,ser,*coeff;
static double logsr2pi=.918938533;
static double b[9]={.035868343,-.193527818,.482199394,-.756704078,
.918206857,-.897056937,.988205891,-.577191652,1.0};

if(x<=0. &&  !((int)((((int)x)-x))) )return errorcode;
/*if( x<0.&& x> -1. )
	{
	return((loggam(1.+x)-log(-x)));
	}
else requires two levels of recursion and  log call,not sin
*/
if (x<-0.) /*was x< -1. when above implemented */
		{/*transform to x>0. will blow up if x integer, as it should*/
		z=1.-x;/* z>2. */
		return(log(pi/ABS(sin(pi*z)))-loggam(z) );
		}
else
	if (x<=1.)/* 0<=x<1 */
		{
		/*z=1.-x*/;/*  0<=z<1*/
		/*return( log(z*pi/sin(pi*z))-loggam(1.+z));*/
		/* Ab& Stegun-takes less than half the time*/
		if(x==0.)return 0.;
		tmp=b[0];
		coeff=&(b[1]);
		for(i=1;i<9;i++)tmp= tmp*x+ *(coeff++);
		return(log(tmp/x));
		}
/* use below for x>1.*/
else
	if(x<=2.)
		{
		tmp=b[0];
		coeff=&(b[1]);
		z=x-1.;
		for(i=1;i<9;i++)tmp= tmp*z+ *(coeff++);
		return(log(tmp));
		}
z=1./x;
tmp=z*z;
/*ser= (1./12.+tmp*(-1./360.+tmp*(1/1260.-tmp/1680.)   ))/x;*/
ser= (.08333333333333+tmp*(tmp*(0.000793650793-.000595238095*tmp)
	-.002777777777))*z;
return (logsr2pi-x+(x-.5)*log(x)+ser);
}

double gamma(x) double x;
{
double y,z;
if(x<=0. &&  !((int)((((int)x)-x))) )return errorcode;
y=exp(loggam(x));
if(x>=0.)return(y);
z=  2*(((int)(-x))%2) -1;
return(y*z);
}

double bernoulli(n) int n;/* n even 0<=n*/
{
double r,sum,term,x,y,b0=1.,b1=-.5,b2=.166666666,g,gamma(double );
int i,itmax=50;
if(n==0)return b0;
if(n==2)return b2;
if(n==1)return b1;
if(n%2==1)return 0.;
r=(double) n;
if(n>3)
	{
	sum=0;
	for(i=1;i<itmax;i++)
		{
		term= pow( (double)(i) ,-r);
		sum+=term;
		if(ABS(term)<5.e-7 || ABS(term/sum)<5.e-7)break;
		}
	g=pow(2.*pi,r);
	x=1.;
	if((n>>1)%2 ==0)x=-1.;
	return 2.*sum*gamma(r+1.)/g*x;
	}
return 0.;
}

