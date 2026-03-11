/* spherical bessel functions and allied routines
C Mathematical Function Handbook
Copyright 1991 Louis Baker. All rights reserved.

sjn(x,n)   returns jn(x) and sets sphb[] array i=0,...n
syn(x,n)   ditto yn(x)
skn(x,n)   ditto kn(x)
sinb(x,n)  all n. (named so as not to conflict with sine[x]=sin(x) )
sj0,sj1,sy0,sy1,si0,si1,sim1,sim2,sk0,sk1- spherical bessel functions
 for n=0,1 (for i, -1 and -2 as well) for j,y,i,k.

for negative n, note: k[x,n]= k[x,-n-1]. e.g. k[x,-1]=k[x,0]
	also y(x,n)= -1^(n+1) j(x,-n-1)
	e. g.  y(x,0) = -j(x,-1)
	Thus j[x,-1] = y[x,0] and y[x,-1]= j[x,0]
*/
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#include "complex.h"
#define max(a,b) ((a)<(b)?(b):(a))
#define min(a,b) ((a)>(b)?(b):(a))

double sj0(x) double x;
{double sin();
if(x==0.)return 1.;
return sin(x)/x;
}
double sy0(x) double x;
{double cos();
if(x==0.)return errorcode;
return -cos(x)/x;
}
double sj1(x) double x;
{double sin(),cos();
if(x==0.)return 0.;
return (sin(x)/x-cos(x))/x;
}
double sy1(x) double x;
{double cos(),sin();
if(x==0.)return errorcode;
return -(cos(x)/x+sin(x))/x;
}

static double y,z,sr;
double si0(x) double x;
{double sinh(),sqrt();
if(x==0.)return 1.;
z=sinh(x)/x;
return z;
}

double sim1(x) double x;
{double cosh(),sqrt();
if(x==0.)return errorcode;
z=cosh(x)/x;
return z;
}
double sim2(x) double x;
{double cosh(),sqrt();
if(x==0.)return errorcode;
z=(sinh(x)-cosh(x)/x)/x;
return z;
}
double si1(x) double x;
{double sinh(),cosh(),sqrt();
if(x==0.)return 0.;
y=(cosh(x)-sinh(x)/x)/x;
return y;
}
double sk0(x) double x;
{double exp(),sqrt(),q;
if(x==0.)return errorcode;
q=pi/(2.*x);
z=q*exp(-x);
return z;
}
double sk1(x) double x;
{double exp(),sqrt(),q;
if(x==0.)return errorcode;
q=pi/(2.*x);
y=q*exp(-x)*(1.+1./x);
return y;
}

double sphb[100];

double sjn(x,n)double x;int n;
{/* as with J, use fwd for x>n else backward*/
int i,m;double old,older,new,sqrt(),s0,s1,mult,syn();
if(n<0)
	{m=-n-1;/* j(x,-1)= -y(x,0)*/
	return (sphb[0]=(m%2?1.:-1.)*syn(x,m));
	}

sphb[0]=s0=sj0(x);
if(!n)return s0;
sphb[1]=s1=sj1(x);
if(n==1)return s1;
if( x>(double)n )
	{
	old=s1;older=s0;sphb[0]=s0;sphb[1]=s1;mult=3.;
	for(i=1;i<n;i++)
		{
		new= old*(mult)/x-older;
		sphb[i+1]=new;mult+=2.;
		older=old;old=new;
		}
	return new;
	}
/* backward*/
if(x==0.)return 0.;
m= min((100-3),n+max(10,sqrt(n*2.)));
/*printf(" m=%d %e\n",m,x);*/
sphb[m+1]=0.;sphb[m]=1.;
for(i=m;i>0;i--)
	sphb[i-1]=((i<<1)+1.)*sphb[i]/x-sphb[i+1];
if(sphb[0]!=0.)new=sj0(x)/sphb[0];
else new=sj1(x)/sphb[1];
/* use when table needed of all n=0 to n */
for(i=0;i<=n;i++)sphb[i]*=new;
return sphb[n];
}

double syn(x,n)double x;int n;
{/*  use fwd */
int i,m;double old,older,new,sqrt(),s0,s1,term,sjn();
if(n<0)
	{m=-n-1;
	return (m%2?1.:-1.)*sjn(x,m);
	}
sphb[0]=s0=sy0(x);
if(!n)return s0;
sphb[1]=s1=sy1(x);
if(n==1)return s1;
	old=s1;older=s0;sphb[0]=s0;sphb[1]=s1;term=3.;
	for(i=1;i<n;i++)
		{
/*		new= old*((i<<1)+1.)/x-older;*/
		new= old*(term)/x-older;
		sphb[i+1]=new;term+=2.;
		older=old;old=new;
		}
	return new;
}

double sinb(x,n)double x;int n;
{/*  use backward */
int i,m,neg,nn;double old,older,new,sqrt(),s0,s1,mult;
sphb[0]=s0=si0(x);
if(!n)return s0;
sphb[1]=s1=si1(x);
if(n==1)return s1;
if(n>0 && x==0.)
	{
	for(i=0;i<=min(n,99);i++)sphb[i]=0.;
	return 0.;
	}

if(n==-1)return sim1(x);
if(n==-2)return sim2(x);
nn=abs(n);
if(n<0)neg=1;
else neg=0;
/* backward*/
m= min((100-3),nn+max(10*(neg+1),sqrt(nn*2.)));
/*printf(" m=%d %e\n",m,x);*/
sphb[m+1]=0.;sphb[m]=1.;
if(neg)
	{
	s1=sim1(x);
	/*if(x<nn)*/
		{/* forward*/
		sphb[0]=s0;sphb[1]=s1;new=1./x;mult=-1.;
		for(i=2;i<=nn;i++)
			{sphb[i]=sphb[i-2]+mult*new*sphb[i-1];mult-=2.;}
		return sphb[nn];
		}
	/*mult= 1.-(m<<1);
	for(i=m;i>0;i--)
		{sphb[i-1]=sphb[i+1]-(mult)*sphb[i]/x;
		mult+=2.;
		} backward recursion not good negative n*/
	}
else
	{/* n>0 backward recursion*/
	mult= (m<<1)+1.;
	for(i=m;i>0;i--)
		{sphb[i-1]=sphb[i+1]-(mult)*sphb[i]/x; mult-=2.;
		}
	}
/*printf(" sphb[0]=%le %le %le\n",sphb[0],sphb[1],sphb[2]);*/
if(sphb[0]!=0.)new=s0/sphb[0];
else new=s1/sphb[1];
/* use when table needed of all n=0 to n */
/*printf(" ck s1=%le new %le,sphb[1] %le rest %le\n",s1,new,sphb[1],new*sphb[1]);*/
for(i=0;i<=nn;i++)sphb[i]*=new;
for(i=1; i<=n;i+=2)sphb[i]=-sphb[i];
return sphb[nn];
}

double skn(x,ni)double x;int ni;
{/*  use fwd */
int i,m,n;double old,older,new,sqrt(),s0,s1,mult;
n=ni;
if(n<0) n= -ni-1;
sphb[0]=s0=sk0(x);
if(!n)return s0;
sphb[1]=s1=sk1(x);
if(n==1)return s1;
	mult=3.;
	old=s1;older=s0;sphb[0]=s0;sphb[1]=s1;
	for(i=1;i<n;i++)
		{
		new=older +old*mult/x;
		mult+=(2.);
		sphb[i+1]=new;
		older=old;old=new;
		}
	return new;
}
