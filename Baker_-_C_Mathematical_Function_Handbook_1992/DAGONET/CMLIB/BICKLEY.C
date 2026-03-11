/*  
Bessel function integrals

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

bii	integrand for adaptive quadrature for bickley functions
bickley	repeated integral of K0
Jrn	repeated integral of Jn

*/
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#define max(a,b) ((a)>(b)?(a):(b))

static double z,r;
extern int level,maxlev,feval,levmax;

double bii(double x)
{double y,ans,num,denom;
if( x==0.) {if(z==0.)
		{if(r==1.)return 2.;
		if(r>1.)return 0.;
		/*r<1.*/
		return pow(2.,r);
		}
	return 0.;}
y= .5*(x+1./x);
num= (z==0.)? 1.: exp(-z*y);
denom= (r==0.)? 1.:pow(y,r);
ans= num/(x* denom );
return ans;
}

/*CAVEAT: should not be used for r <1. Strictly speaking,
r should be an integer>= 1*/

double bickley(double x, double rr)
{
double a,b,adsimp(),ans,dx,cutoff;int i,lim, lm;
r=rr;z=x; maxlev= 9/*6*/;cutoff=.00015625*.5;
lim=3;if(rr<2.)lim=8/max(rr,.01);ans=0.;
a=0.;dx=cutoff/lim;b=dx;lm=0;
for(i=0;i<lim;i++,a+=dx,b+=dx)
	{
	ans+= adsimp(a,b,1.e-7,bii);
	lm=max(lm,levmax);
	}
lim=3;dx= (1.-cutoff)/lim;a=cutoff;b=a+dx;
for(i=0;i<lim;i++,a+=dx,b+=dx)
	{
	ans+= adsimp(a,b,1.e-7,bii);
	lm=max(lm,levmax);
	}
	if( lm>=maxlev)
		fprintf(stderr,
			" warning: bickley may not have achieved desired accuracy\n");

return ans;
}


#define tol 1.e-7

extern double bessela[100];

#define FILL 25

double Jrn(double x, double r, int n)
{int method,k,intnr,id,lim=99;double sum,term,factor,jj,rv;
/*printf(" entered Jrn\n");printf(" x=%le r=%le n=%d\n",x,r,n);*/
method=r; if(abs(r-method)<tol && r<=0.)return errorcode;
method=1;sum=n+r;intnr=(int)sum;
if( (sum-(double)intnr) ==0. && intnr>=0)
	{method=0;
	sum=jn(x,FILL); /*set up bessela array*/
	/*for(k=0;k<100;k++)
		{printf(" J[%d]=%le\n",k,bessela[k]);}*/
	}

/*printf(" method=%d\n",method);*/

for(factor=1.,sum=0.,k=0;k<lim;k++)
	{

	if(method)jj=jbes(x,r+(double)(n+(k<<1)));
	else { id=intnr+(k<<1);if(id> FILL-1 )break; jj=bessela[id];}
	term=gamma((double)k+r)*factor*jj;
	sum+=term;
	if(abs(term)<abs(sum)*tol)break;
	factor/= (k+1.);
	}
rv=pow(2.,r)/gamma(r)*sum;
return rv;
}
