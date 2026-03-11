/* Mehler (Conical Legendre) functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

Mehler(x,z,mm) Mehler conical function P  [ix-.5] (z)
Mehler0(x,z) double x,z;
*/

#include "cmlib.h"
#include "protom.h"
#include <stdio.h>
#define iabs(x) ((x)<0?-(x):(x))
#define tol 1.e-9

double Mehler(x,z,mm) double x,z;int mm;
{/* P ix-1/2(z) */
double y,sum,t,q,term,factor,pow(),kmm,gamma(),tan,kk,power;
int m,k;
/* z=cos theta, y =sin(theta/2)^2*/
if(z==-1.)return errorcode;
y=.5*(1.-z);tan= (1.-z)/(1.+z);
/* k -m,x(z)*/ m=iabs(mm);
sum=1.;t= 4.*x*x;q=1.;term=t+1.;
factor=term/(4.*(1.+m))*y;
for(k=1;k<1000;k++)
	{sum+=factor;
	q+=2.;
	term=t+q*q;
	factor*=term/((k+1.)*4.*(m+1.+k))*y;
	if(abs(factor)<tol*abs(sum))goto fini;
/*printf(" Mehler sum=%e latest term %e\n",sum,factor);*/
	}
fprintf(stderr," Mehler tolerance not achieved\n");
fini:power= m?pow(tan,(double)m):1.;
kmm= (m%2?-1.:1.)/gamma(m+1.)*power*sum;
if(mm<=0)return kmm;
kk=kmm;t*=.25;/*t=x*x*/
for(k=1;k<=(m<<1)-1;k+=2)
	{
	kk*=(t+.25*(k*k));
	}
return kk;
}

double Mehler0(x,z) double x,z;
{/* P ix-1/2(z) */
double y,sum,t,q,term,factor,pow(),gamma(),tan,cosh(),log(),g,yy,lt,r,ss;
int k;
/* z=cos theta, y =sin(theta/2)^2*/
y=.5*(1.-z);t=4.*x*x;term=t+1.;q=1.;
/* k m=0,x(z)*/
if(z==-1.)return errorcode;
if(y>.95)
	{/* near log. singularity at theta=pi, cos()=-1,y=1*/
	yy= .5*(1.+z);tan= (1.-z)/(1.+z);lt=log(tan);
	for(g=0.,k=1;k<100;k++)
		{
		q=(k<<1)-1;
		term= (t-q)/(k*(t+q*q));
		g+=term;
		if(abs(term)<tol*abs(g))goto fini;
		}
		fprintf(stderr," Mehler0(y>.95 1st) tolerance not achieved\n");
		fini: g*=2.;
/*printf(" g= %e cos^2 theta/2=%e\n",g,yy);*/
	factor=(t+1.)/(4.*(1.))*yy;sum=0.;
	ss=0.;r=1.;q=1.;
	for(k=1;k<1000;k++)
		{
		ss+=1./r;r+=1.;
		/*printf(" ss=%e\n",ss);*/
		tan=factor*(lt-g+2.*ss);
		sum+=tan;
		/*printf(" sum=%e factor %e latest term %e\n",sum,factor,tan);*/
		if(abs(tan)<tol*abs(sum))goto fin;
		q+=2.;
		term=t+q*q;
		factor*=term*yy/((k+1.)*4.*(1.+k));
		}
		fprintf(stderr," Mehler0(y>.95 2nd) tolerance not achieved\n");
		fin: return cosh(pi*x)/pi*(lt-g+sum);
	}
sum=1.;t= 4.*x*x;q=1.;term=t+1.;
factor=term/4.*y;
for(k=1;k<1000;k++)
	{sum+=factor;
	if(abs(factor)<tol*abs(sum))goto fin2;
	q+=2.;
	term=t+q*q;
	factor*=term/((k+1.)*4.*(1.+k))*y;
	}
fprintf(stderr," Mehler0 tolerance not achieved\n");
fin2: return sum;
}
