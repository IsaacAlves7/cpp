/*  
iterated error function. backward recursion used, as
well as forward for small arguments.  

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <alloc.h>
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"


static double *tierfc=NULL;

double ierfc(z,n) int n; double z;
{/* uses stable backward recursion for large z. this give poor results for
small z, however*/
int m,i;
double *w,norm,coef,ans,ratio;
if(!n){ norm=erf(z,&coef); if(tierfc)*tierfc=coef; return coef;}/* erfc*/
if(n==-1) return exp(-z*z)*2./sqrt(pi); 
if(z<0. || n<-1)return errorcode;
if(z<.9)return ierfcf(z,n);
m=n+15;
if(z<2.)m+=10;
if(z<1.75) m+=10;
if(z<1.5) m+=10;
if(z<1.)m+=10;
w=(double *)malloc( (m+3)*sizeof(double));
if(w==NULL)
	{
	fprintf(stderr," ierfc could not allocate space for array");
	return errorcode;/* or 0., as the result is quite small*/
	}
w[m+1]=1.;w[m+2]=0.;
for(i=m;i>=-1;i--)
	{
	norm= 2.*(i+2)*w[i+2]+z*2.*w[i+1];
	/*printf(" w[%d]= %e\n",i,norm);*/
	if (i!=-1) w[i]=norm;
	}
coef=2./sqrt(pi)*exp(-z*z);
/*printf(" coef=%e\n",coef);*/
ratio=coef/norm;
ans=w[n]*ratio;
if(tierfc)
	{for(i=0;i<=n;i++) tierfc[i]=w[i]*ratio;}
free(w);
return ans;
}


double ierfcf(x,n) double x;int n;
{/* uses forward recursion-can be unstable. not recommended*/
double erf(),y,z,exch;
int k;
y=erf(x,&z);
if(tierfc)*tierfc=z;
if(!n)return z;
y=2./sqrt(pi)*exp(-x*x);
if(n==-1)return y;
if(z<0. || n<-1)return errorcode;
/* y=i^-1,z=i^0*/
for(k=1;1;k++)
	{
	y=(.5*y-x*z)/k;
	if(tierfc)tierfc[k]=y;
	if(k==n)break;
	exch=z;
	z=y;
	y=exch;
	}
return y;
}


/* returns table[k]= i ^k erfc (z) for k=0,...,n*/
ierfctable(z,n,table) int n; double z, table[];
{double dum;
tierfc=table;
dum=ierfc(z,n);
tierfc=NULL;
if(dum==errorcode){table[0]=errorcode;return ierrorcode;}
return 0;
}
