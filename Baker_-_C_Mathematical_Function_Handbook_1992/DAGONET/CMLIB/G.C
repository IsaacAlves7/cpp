
/*  
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

loggam,gamma	logarithm of gamma, gamma functions all real x.
P incgam BigGamma lcgam, ucgam
		versions of incomplete gamma function
fac		factorial
pochhammer	pochhammer symbol for real x, integer n
cgamma		gamma, log gamma for complex arguments
gammaqd		compact, "exact" gamma.
cdigamma	digamma for complex arguments
polygam		polygamma function for real arguments
*/

#include <stdio.h>
#include "cmlib.h"
#define DOFOR(I,J) for(I=0;I<J;I++)
#include "complex.h"
#include "protom.h"

#define min(a,b) (((a)<(b))? (a): (b))

int itkt;
int iterp;/* global to return count used*/


double gamma(x) double x;
{
double y,z;
if(x<=0. &&  (((x-(double)((int)x)))==0.) )
	{
	fprintf(stderr," gamma: arg 0 or neg. integer=%le\n",x);
	return errorcode;
	}
y=exp(loggam(x));
if(x>=0.)return(y);
z=  2*(((int)(-x))%2) -1;
return(y*z);
}


double loggam(x) double x;
{
int i;
double z,tmp,ser,log(),sin(),*coeff;
static double logsr2pi=.918938533;
static double b[9]={.035868343,-.193527818,.482199394,-.756704078,
.918206857,-.897056937,.988205891,-.577191652,1.0};

/*if( x<0.&& x> -1. )
	{
	return((loggam(1.+x)-log(-x)));
	}
else requires two levels of recursion and  log call,not sin
*/
if (x<-0.) /*was x< -1. when above implemented */
		{/*transform to x>0. will blow up if x integer, as it should*/
		z=1.-x;/* z>2. */
		return(log(pi/abs(sin(pi*z)))-loggam(z) );
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
ser= (1./12.+tmp*(-1./360.+tmp*(1/1260.-tmp/1680.)   ))*z;
/*ser= (.08333333333333+tmp*(tmp*(0.000793650793-.000595238095*tmp)
	-.002777777777))*z;*/
return (logsr2pi-x+(x-.5)*log(x)+ser);
}
#define small 1.e-30

double P(a,x) double x,a;
{/* P incomplete gamma function*/
int i,itmax=100;
double gln,exp(),log(),loggam(),sum,ap,del,fi,start,
tol=3.e-7,c0,d0,ana,offset,mult,delta;
/* error condition return -1 on invalid arguments*/
if( x< 0. || a<0. ) return(-1.);
if(x==0.)return(0.);
gln=loggam(a);
if (x< (a+1.))
	{
	/*series*/
	offset=0.;
	mult=1.;
	ap=a;
	sum=1./a;
	del=sum;
	DOFOR(i,itmax)
		{
		ap++;
		del*=x/ap;
		sum+=del;
		if( abs(del)<abs(sum)*tol) goto fini;
		}
		fprintf(stderr," trouble incomplete gamma series\n");
	}
else
	{
	offset=1.;
	mult=-1.;
	start=small;
	sum=start;/*fictitious should be small but numerical pblms?*/
	d0=0.;c0=sum;
	DOFOR(i,itmax)
		{
		fi=i;
		if(i)ana=fi;
		else ana=1.;
		d0= (x+d0*ana);
		c0=(x+ana/c0);
		if(d0==0.)d0=small;
		if(c0==0.)c0=small;
		d0=1./d0;
		delta=d0*c0;sum*=delta;
		ana=fi+1.-a;
		d0= (1.+d0*ana);
		c0=(1.+ana/c0);
		if(d0==0.)d0=small;
		if(c0==0.)c0=small;
		d0=1./d0;
		delta=d0*c0;sum*=delta;
		if( abs(delta-1.)<tol)
			{/*sum-=start;*/ goto fini;}
		}
	fprintf(stderr," trouble incomplete gamma cont. fract\n");
	}
/*return(-1.);*/
fini:return(offset+mult*sum*exp(-x+a*log(x)-gln));
}
double incgam(a,x) double a,x;
{return P(a,x);}
double lcgam(a,x) double a,x;
{double gamma();return P(a,x)*gamma(a);}
double ucgam(a,x)double a,x;
{double gamma(); return gamma(a)*(1.-P(a,x));}

double BigGamma(a,x) double x,a;
{/* incomplete gamma function*/
int i,itmax=100;
double exp(),log(),gamma(),mult,offset,sum,ap,del,fi,start,
tol=3.e-7,c0,d0,ana,delta;
if( x< 0. || a<0. ) return(errorcode);
if(x==0.)return(0.);
if (x< (a+1.))
	{
	/*series */
	offset=gamma(a);
	mult=-1.;
	ap=a;
	sum=1./a;
	del=sum;
	DOFOR(i,itmax)
		{
		ap++;
		del*=x/ap;
		sum+=del;
		if( abs(del)<abs(sum)*tol) goto fini;
		}
		fprintf(stderr," trouble incomplete gamma series\n");
	}
else
	{
	offset=0.;
	mult=1.;
	start=small;
	sum=start;/*fictitious should be small but numerical pblms?*/
	d0=0.;c0=sum;
	DOFOR(i,itmax)
		{
		fi=i;
		if(i)ana=fi;
		else ana=1.;
		d0= (x+d0*ana);
		c0=(x+ana/c0);
		if(d0==0.)d0=small;
		if(c0==0.)c0=small;
		d0=1./d0;
		delta=d0*c0;sum*=delta;
		ana=fi+1.-a;
		d0= (1.+d0*ana);
		c0=(1.+ana/c0);
		if(d0==0.)d0=small;
		if(c0==0.)c0=small;
		d0=1./d0;
		delta=d0*c0;sum*=delta;
		if( abs(delta-1.)<tol)
			{/*sum-=start;*/ goto fini;}
		}
	fprintf(stderr," trouble incomplete gamma cont. fract\n");
	}
fini:return(offset+mult*sum*exp(-x+a*log(x)));
}

double Pgamma(a,x)double a,x;
{double gamma(); return BigGamma(a,x)/gamma(a);}

double Smallgamma(a,x)double a,x;
{return P(a,x)*gamma(a);}

double fac(x) double x;
{
double gamma();
return gamma(x+1.);
}

double pochhammer(z,n) int n;double z;
{double gamma(),m;int k;
if(!n)return 1.;
k=z;
if( (z-k)!=0. || z>0.)
	return gamma(z+n)/gamma(z);
/* z is a negative integer*/
if(n==1)return z;
m=z;
for (k=1;k<n;k++)
	{m *=(z+k);if(m==0.)return 0.;}
return m;
}

/* quick & dirty gamma function*/
#define lim 20

double gammaqd( double xin)
{double xi,factor=1.,g,x,sin(),exp(),sqrt(),pow();
/*reflect*/
x=xin;
if(x==0.)return errorcode;
if(x<0.)
	{
	if(abs( x-(int)x)<1.e-10)return errorcode;
	return pi/(gamma(1.-x)*sin(pi*x));
	}
/* may take x>0 at this point*/
while(x< lim){factor=factor/x;x++;}xi=1./x;
g= factor*(exp(-x))*(pow(x,x-.5))*sqrt(2.*pi)*
(1.+xi*(1./(12.)+xi*(1./(288.)-xi*(139./(51840.)+571./(2488320.*x)))));
return g;
}

/* complex gamma*/


cgamma(zz,ans,loggam) struct complex *zz,*ans,*loggam;
{struct complex z,zm,t,tt,sum,term,den,a,aux;
 double c[12]={1./12.,-1./360.,1./1260.,
-1./1680., 1./1188.,-691./360360.,1./156.,
 -3617./122400.,43867./244188.,-174611./1125400., 77683./5796.},
  x,y,tol=1.e-7,xdist,log();
 int flip;
 /* best of the bunch- from CALGO 404 Lucas and Terrill */
 /* modified 9/25/90 so that |gamma(x+iy)|=|gamma(x-iy)|
 independently of how argmt() returns answer*/
 int i,m,reflect;
 CLET(z,*zz);flip=0.;
 x=z.x;y=z.y; reflect=0;
 if( x<tol)
	{
	xdist=x- ((int)(x-.5));
	CMPLX(zm,xdist,y);
	if(cabs(zm)<tol)
		{ans->x=errorcode; ans->y=0.;return 0;}
	if(x<0.)
		{
		reflect=1;x=1.-x;y=-y;CMPLX(z,1.,0.);CSUB(z,z,*zz);
		}
	}
 if(y<0.)
	{flip=1;y=-y;z.y=y;}
 m=0; 
  while(x<10.)
	{
	x+=1.;
	m++;
	}
 while( /*abs*/(y)>=x)
	{
	x+=1.;m++;
	}
 CMPLX(t,x,y);
 CMULT(tt,t,t);CLET(den,t);
 CMPLX(aux,.5,0.);CSUB(aux,t,aux);
 clog(&t,&zm);CMULT(sum,aux,zm);CSUB(sum,sum,t);
 sum.x+= .5*log(2.*pi);
 for(i=0;i<11;i++)
	{
	CMPLX(aux,c[i],0.);CDIV(term,aux,den);
	if(abs(term.x)<abs(sum.x)*tol)
		{
		if(y==0. ||(abs(term.y)>=abs(sum.y) ) )break;
		}
	CADD(sum,sum,term);
	CMULT(aux,den,tt);CLET(den,aux);
	}
 if(m)
	{for(i=0;i<m;i++)
		{CMPLX(a, (double)i,0.);CADD(a,a,z);
		clog(&a,&aux); CSUB(sum,sum,aux);
		}
	}
 if(reflect)
	{
	CMPLX(aux,pi,0.);
	CTREAL(t,z,pi);
	csin(&t,&tt);
	CDIV(den,aux,tt);
	clog(&den,&aux);
	CSUB(sum,aux,sum);
	}
 if(flip) {sum.y=-sum.y;}
 CSET(loggam,sum);cexp(&sum,ans);
/*printc(loggam);printc(ans);printf(" before flip\n");*/
/* if(flip) {ans->y=-ans->y;loggam->y=-loggam->y;}*/
 /* flip sign of imaginary part of gamma. as log(x)=|z|+i arg(z),
 flip sign of imag. part of logarithm*/
 return 0;
 }



extern double digammin;

cdigamma(x,ans) struct complex *ans,*x;
{struct complex q,y,z,sum,rq,one;
if(x->x==1. && x->y==0.){CMPLX(*ans, - Egamma ,0.);return 0;;}
CLET(q,*x);CMPLX(sum,0.,0.);CMPLX(one,1.,0.);
if(q.x < 0.){
	if( abs((int)(q.x)-q.x && q.y==0.) <1.e-8)
		{CMPLX((*ans),errorcode,errorcode);return 1;}
	CSUB(q,one,*x);
	CTREAL(rq,*x,pi);
	ccot(&rq,&sum);
	CTREAL(sum,sum,pi);
	}
	/*note sum subtracted from final result*/
if(digammin<20.)digammin=20.;
while(q.x<digammin)
	{
	if(q.x==0. && q.y==0.){CMPLX((*ans),errorcode,errorcode);return 1;}
	CDIV(rq,one,q);CADD(sum,sum,rq);
	q.x+=1.;
	}
CDIV(y,one,q);/*y=1/q*/;CMULT(z,y,y);clog(&q,&rq);
/* logx-1/2x- sum n=1 to inf of  B[2n]z^-2n/2n*/
CSUB(*ans,rq,sum);CTREAL(rq,y,-.5);CADD(*ans,*ans,rq);
CLET(sum,(z));CTREAL(sum,sum,-1./240.);
sum.x+=1./252.;CMULT(rq,sum,z);rq.x-= 1./120.;
CMULT(sum,rq,z); sum.x+=1./12.;	CMULT(rq,sum,z);
CSUB(*ans,*ans,rq);
return 0;
}

#define itmax 1000

double polygam(z,n) double z; int n;
{/* wrong for n=0. better to use digam, etc. where possible*/
double x; int m; m=n+1; x=(double)m;
return ( (n%2)? 1.:-1.) *gamma(x)*zeta2( x,z);
}
