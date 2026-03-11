/*  
integrals of bessel functions 
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

ji	integral of J0
jii	1. - integral of J0
jin	integral of Jn
jiotn	integral of Jn/t dt
ki	integral of K0 from 0 to x
kii	pi/2 - integral of K0 above = integral of K0 from x to infinity
ij0t	integral (1-J0_/t dt
ii0m1t	integral (I0 - 1)/t
iy0t	integral Y0/t dt
ik0t	integral K0/t dt
y0i	integral of y0	uses Struve functions
i0i	integral of i0   "     "        "

*/

#include "cmlib.h"
#include "protom.h"

#define gamma   Egamma
#define eps	3.e-7
#define epsk	1.e-8

/*static int itkt;*/
extern int itkt;

double ji(x) double x;
{
double y,p,q,r,ans;
y=j0(x);
p=j1(x);
q=h0(x);
r=h1(x);
return x*(y+.5*pi*(q*p-r*y));
/*return x*(y+.5*pi*(h0(x)*j1(x)-h1(x)*y));*/
}

double jii(x) double x;
{
double ji(),z;
return 1.-ji(x);
}

double jin(x,n) int n;double x;
{/* integral from 0 to x of j[n](x)*/
int k,m;
double sum,jn(),ji();
sum=0.;
m=n>>1;
if( n%2)
	{/* n odd*/
	if(n==1) return 1.-j0(x);
	for(k=1;k<=m;k++)
		sum+=jn(x,k<<1);
	return 1.-j0(x)-2.*sum;
	}
for(k=0;k<m;k++)
	{
	sum+= jn(x, (k<<1)+1 );
	}
return ji(x)-sum*2.;
}

double jiotn(x,n) int n;double x;
{
int k,m,halfn;
double sum,jn(),y,z,ji(),j1();
if(n<1)return errorcode;
sum=0.;
halfn=n>>1;
if( n%2)
	{/* n odd*/
	if(x==0.)return 0.;
	if(n==1) return ji(x)-j1(x);
	for(k=1;k<=halfn;k++)
		sum+=jn(x,k<<1)*k;
	return (ji(x)-j1(x)-4./x*sum)/(2*halfn+1);
	}
for(k=1;k<=halfn;k++)
	{
	m=(k<<1)-1;
	sum+= m*jn(x, m );
	}
return (1.-sum*2./x)/(2*halfn);
}

double ki(x) double x;
{/* integral from 0 to x of K0*/
double y,exp(),sqrt(),l0(),k0(),l1(),k1(),p,q,r,ans;
if(x==0.)return 0.;
if(x>7.)
	{
	 y=7./x;
	 return pi*.5-exp(-x)/sqrt(x)*
	 ((((((.00033934*y-.00163271)*y+.00417454)*y
	 -.00933994)*y+.02576746)*y-.11190289)*y+1.25331414);
	}
/* x<7*/
/* a triple sum! A&S p 480*/
/* or use  11.1.8 with K0,K1,L0,L1 L the Modifed Struve functions*/
y=k0(x);
p=k1(x);
q=l0(x);
r=l1(x);
ans= x*(y+pi*.5*(y*r+q*p ));
return ans;
}


double kii(x) double x;
{
double ki(),z;
return pi*.5-ki(x);
}



double ii0m1t(x) double x;
{
int k;double sum,term,factor,y,fk,ii0tas();
if(x>5.)return ii0tas(x);
sum=0.;
term=1.;factor=1.;
y=x*x;
for(k=1;k<50;k++)
	{
	fk=k;
	factor*= y*.25/(fk*fk);
	term=factor/k;
	sum+=term;
	if(abs(term)<eps || abs(sum)*eps>abs(term))break;
	}
itkt=k;
return sum*.5;
}

double ij0t(x) double x;
{      int k;
double log(),y,z,sum,term,fk,factor,digamma();
y=.5*x;
sum=0.;y*=y;factor=1.;fk=1.;
for(k=1;k<50;k++)
	{
	factor*=-y/(fk*fk);
	term=factor/fk;
	sum+=term;
	if( abs(sum)*eps>abs(term))break;
	fk+=1.;
	}
itkt=k;
return -.5*sum;
}

double iy0t(x) double x;
{      int k;
double log(),y,z,sum,term,fk,factor,digamma();
y=.5*x;
z=log(y);
sum=0.;y*=y;factor=1.;fk=1.;
for(k=1;k<50;k++)
	{
	factor*=-y/(fk*fk);
	term=factor/fk*(digamma(1.+fk)+.5/fk-z);
	sum+=term;
	if( abs(sum)*eps>abs(term))break;
	fk+=1.;
	}
itkt=k;
return -z*z/pi-2./pi*gamma*z+(pi*pi/6.-gamma*gamma+sum)/pi;
}

double ik0t(x) double x;
{      int k;
double log(),y,z,sum,term,fk,factor,digamma(),ik0tas(),digam();
/*if(x>5.)return ik0tas(x);*/
y=.5*x;
z=log(y);
sum=0.;y*=y;factor=1.;fk=1.;
for(k=1;k<50;k++)
	{
	factor*=y/(fk*fk);
	term=factor/fk*(digamma(1.+fk)+.5/fk-z);
	sum+=term;
	if(abs(sum)*epsk>abs(term))break;
	fk+=1.;
	}
itkt=k;
return z*z*.5+gamma*z+pi*pi/24.+gamma*gamma*.5-.5*sum;
}

double ii0tas(x) double x;
{
double sqrt(),exp(),y,sum;int i;
double b[11]={
.3989314,.1332055,-.0493843,1.4780044,-8.6556013,28.1221478,-48.0524115,
40.3947340,-11.9094395,-3.5195009,2.1945464};
if(x<5.)return errorcode;sum=b[10];y=5./x;
for(i=9;i>=0;i--) sum=(b[i]+sum*y);
return  sum/(sqrt(x)*x*exp(-x));
}

double ik0tas(x) double x;
{
double sqrt(),exp(),y,sum;int i;
double b[7]={1.2533141,-.5091339,.3219184,-.2621446,.2060126,-.1110396,.02724};
if(x<4.)return errorcode;
sum=b[6];y=4./x;
for(i=5;i>=0;i--) sum=(b[i]+sum*y);
return  sum/(sqrt(x)*x*exp(x));
}


double y0i(x) double x;
{double h0(),h1(),y0(),y1(),y;
y=y0(x);
return x*(y+pi*.5*(h0(x)*y1(x)-h1(x)*y));
}

double i0i(x) double x;
{double i,ii,s0,s1;struct complex ans,arg;
/*i=i0(x);
return x*(i+pi*.5*(l1(x)*i-l0(x)*i1(x)));*/
/*i=in(x,0);
return x*(i+pi*.5*(l1(x)*i-l0(x)*in(x,1)));*/
CMPLX(arg,x,0.);
ibess(&arg,&ans,0);
i=ans.x;
ibess(&arg,&ans,1);
ii=ans.x;
s0=StruveL(0.,x);
s1=StruveL(1.,x);
/*printf(" i0 %le i1 %le L0 %le L1 %le\n",i,ii,s0,s1);*/
return x*(i+pi*.5*(s1*i-s0*ii));
}
