/* digamma and first 2 polygamma functions
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

digam	=(d/dx) (log Gamma). 
	former rational approx.,latter 'exact'

digamma,pg1,pg2	1st,2nd deriv. of digamma
	based on approx. in Tuma. |error|< 1.e-10

*/

#include "cmlib.h"
#include "protom.h"

double digammin;

double digamma(x) double x;
{double q,y,z,log(),sum,tan();
q=x;    sum=0.;
if(q==0.)return errorcode;
if(digammin<=0.)digammin=30.;
if(q<0.){
	if( abs((int)(q)-q) <1.e-8)return errorcode;
	q=1.-x;sum=pi/tan(pi*x);}/*note sum subtracted from final result*/
while(q<digammin)
	{
	if(q==0.)return errorcode;
	sum+=1./q;
	q++;
	}
y=1/q;z=y*y;
/* logx-1/2x- sum n=1 to inf of  B[2n]z^-2n/2n*/
return log(q)-sum-.5*y-z*(1./12.-z*(1/120.-z*(1/252.-z/240.)));
}

double pg1(x) double x;
{
int i;
double sum,z,u,sin(),y;
double b[10]={.4041138064,.1471100206,.0500956643,
.0160671426,.004941886,.0014725602,.0004282353,
.0001221952,.000034347,.0000095387};
if( x<=0. && abs((int)(x)-x) <1.e-8)return errorcode;
u=x-1.;
if(u>0. && u<=.5)
	{
	sum=b[9];
	y=u*u;
	for(i=8;i>=0;i--)	
		{
		sum= sum*y+b[i];
		}
	sum*=u;
	/*printf(" omega=%e\n",sum);*/
	z=1./(u+1.); /* z=1/x */
	sum-=.5*z*z;
	z=1./(1.-u);
	sum+=.5*z*z;
	z=sin(u*pi);
	return  pi*pi*.5/(z*z)-sum-.5/(u*u);
	}
if(x>=2.)
	{
	/*printf(" asymptotic x=%e\n",x);*/
	y=1./x;
	z=y*y;
/*printf(" asymptotic z,y=%e %e\n",y,z);*/
/* desmet blows up on underflow*/
	return(y+z*(.5+y*(1./6.+z*(-1./30./*+z*(1./42.-z/30.)*/))));
	}
if(x<0.)
	{
	z=sin(pi*x);
	return (-pg1(1.-x)-pi*pi/(z*z));
	}
/*  .5<=x<2.*/
/*if(x>1.)*/
return(pg1(x+1.)+1./(x*x));
}

double pg2(x)double x;
{
int i;
double sum,z,u,y,sin(),cos();
double b[10]={.4041138064,.4431326652,
.250478322,.1124699968,.04476974,
.0161981556,.0055670524,.001832928,.0005839092,.000181222};
if( x<=0. && abs((int)(x)-x) <1.e-8)return errorcode;
u=x-1.;
if(u>0. && u<=.5)
	{
	sum=b[9];
	y=u*u;
	for(i=8;i>=0;i--)
		{
		sum= sum*y+b[i];
		}
	z=1./(u+1.);
	sum+=z*z*z;
	z=1./(1.-u);
	sum+=z*z*z;
	z=sin(u*pi);
	return  -pi*pi*pi*cos(pi*u)/(z*z*z)-sum+1./(u*u*u);
	}
if(x>=2.)
	{
	y=1./x;
	z=y*y;
	/* prevent underflows for DeSmet*/
	return(-z+z*(-y+z*(.5+z*(1./6./*+z*(-1./6.+z*(.3-5.*z/6.))*/))));
	}
if(x<0.)
	{
	z=sin(pi*x);
	return (pg2(1.-x)+pi*pi*pi*cos(pi*z)/(z*z*z));
	}
/*  .5<=x<2.*/
/*if(x>1.)*/
return(pg2(x+1.)-2./(x*x*x));
}


double digam(x) double x;
{
/* caveat- Tuma defines ad dlog gamma(x+1)/dx,
A & S d log gamma(x)/dx
*   will conform to A&S*/
double sum,z,y,g,log(),tan();
double b[9]={ Egamma ,.2020569031,.0369277551,
.0083492774,.0020083928,.0004941886,.0001227133,
.0000305882,.0000076372};
double c[11]={0.,0.,.4227843351,.9227843351,1.2561176684,1.5061176684,
1.7061176684,1.8727843351,2.0156417780,2.1406414780,2.25175258918};
double u;int i;
if( x<=0. && abs((int)(x)-x) <1.e-8)return errorcode;
g= Egamma;
if(x==1.)return -g;
if(x<=10.)
	{i=x; if( abs(x-(double)i)<1.e-8)return c[i];	}
u=x-1.;
if(u>-0.5 && u<=.5)/*x between .5 and 1.5*/
	{
	sum=b[8];
	y=u*u;
	for(i=7;i>=0;i--)
		{
		sum= sum*y+b[i];
		}
	if(u>0.)
	return 1.+.5/u- 1./((1.+u)*(1.-u))-pi*.5/tan(pi*u)-sum;
	return 1.-.5/u- 1./((1.+u)*(1.-u))+pi*.5/tan(pi*u)-sum;
	}
if(x>2.)
	{
	y=1./x;
	z=y*y;
	return(log(x)-y*.5+z*(z*(1./120.-z/252.)-1./12.));
	}
if(x<-.5)
	return (digam(1.-x)-pi/tan(pi*x));
/*  .5<=x<2.*/
/*if(x>1.)*/
return(digam(x+1.)-1./(x));
/* .0<=x<=.5 also*/
}
