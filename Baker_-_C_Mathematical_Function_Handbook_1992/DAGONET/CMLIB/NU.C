/* nu and mu integrals of Erdelyi Vol III p.217
mu(x,alpha,beta) --caveat: Erdelyi has mu(x,beta,alpha)
nu=mu(x,alpha,0) or mu(x,0,0)

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include "cmlib.h"
#include "protom.h"

#define max(a,b) ((a)>(b)?(a):(b))

static double z,co,alphan,betan;
extern int level,maxlev,feval,levmax;
int lim1,lim2,lim3;

double nui(t) double t;
{double ans;

ans=pow(z,t)/gamma(1.+t+alphan);
if(betan!=0.)ans*=pow(t,betan);
return ans;
}


#define toll 1.e-10

double mu(double x,double ain,double bin)
{
double a,b,ans,dx,cutoff,print;int i,lim;
if(x<=0.)return errorcode;
/*printf(" in mu %le %le %le %d\n",x,ain,bin,lim1);*/
co=-1.;lim1=max(lim1,1);
z=x; maxlev= 6;cutoff=co;
if(co<=0.)
	{/*automatic selection of co. note bi>0 x>0.*/
	a=max(1.,nui(0.));cutoff=4.;
	infinite_loop
		{
		dx=nui(cutoff);
		if( dx< toll*a)break;
		cutoff*=2.;
		}
	}
co=cutoff;
/*printf(" cutoff=%le\n",cutoff);*/
a=0.;ans=0.; alphan=ain;betan=bin;
lim=lim1;dx= cutoff/lim;b=a+dx;
for(i=0;i<lim;i++,a+=dx,b+=dx)
	{
	/*printf(" a,b %le %le ans=%le now\n",a,b,ans);*/
	ans+=/* print=*/adsimp(a,b,1.e-7,nui);
	/*printf(" a,b,dans=%le %le %le now\n",a,b,print);*/
	}
/*printf(" ans after 1- %le\n",ans);*/
if(ain!=0.)ans*=pow(x,ain);
if(bin!=0.)ans/=gamma(1.+bin);
return ans;
}

double nu(double x)
{return mu(x,0.,0.);}

double Nu(double x,double a)
{return mu(x,a,0.);}

