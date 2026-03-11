
/* "Abramowitz" functions f1, f2, f3 from Ch.27 

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
extern int itkt;


double af2s(x) double x;
{
double log(),ln,sqrt(),sp,sq,sum1,sum2;
sp=sqrt(pi);
if(x>0.)ln=log(x);
else ln=0.;
sq=x*x;
sum2=.3333333+sq*(.000198*sq-.01667);
sum1= sp*.5*(1.+sq)+x*(sq*(-.3225+x*(-.1477
+x*(.03195+x*(.00328+x*(-.000491-.0000235*x)))))-1.);
/*sum1= .5*sp-x+.5*sp*sq-.3225*cube-.1477*sq*sq
+.03195*cube*sq+.00328*cube*cube-.000491*cube*cube*x
-.0000235*cube*cube*sq;*/

return .5*(ln*x*sq*sum2+sum1);
}

double af3s(x) double x;
{
double log(),ln,sqrt(),sp,sq,sum1,sum2;
sp=sqrt(pi);
if(x>0.)ln=log(x);
else ln=0.;
sq=x*x;
sum2=.0833+sq*(-.00278+.000025*sq);
/*sum1= 1.-.5*sp*x+.5*sq-.2954*cube+.1014*sq*sq
+.02954*cube*sq-.00578*cube*cube-.00047*cube*cube*x
+.000064*four*four;*/
sum1= 1.+x*(-.5*sp+x*(.5+x*(-.2954+x*(.1014+x*(
.02954+x*(-.00578+x*(-.00047+.000064*x)))))));

return .5*(sum1-ln*sq*sq*sum2);
}

double af1s(x) double x;
{
int k,itmax=400,in;
double log(),sqrt(),sp,sq,ln,sum,term,coef,a,b,aa[2],bb[2],t;
sp=sqrt(pi);
sq=x*x;
if(x>0.)ln=log(x);
else ln=0.;
sum=0.;
term=1.;
in=0;
aa[0]=0.;aa[1]=-1.;
bb[0]=-sp;bb[1]=1.5*(1.-.57721556649);
a=0.;b=1.;
for(k=0;k<itmax;k++)
	{
	coef= a*ln+b;
	t=coef*term;
	sum+=t;
/*printf(" k in %d %d a,b %e %e s t %e %e c=%e\n",k,in,a,b,sum,t,coef);*/
	if( abs(t)<1.e-6 || abs(t)<1.e-5*abs(sum))break;
	term*=x;
	itkt=k;
	if(k>1)
		{
		coef=1./((k+1)*(k-1)*k);
		aa[in]= -2.*aa[in]*coef ;
		bb[in]= (-2.*bb[in]-(3*k*k-1)*aa[in])*coef;
		}
	a=aa[in];b=bb[in];
	in++;
	if(in>1)in=0;
	}
return .5*sum;
}


double afa(x,n) double x,n;
{
int i,itmax=30,k;
double sqrt(),exp(),pow(),y,z,sum,term,v,vi,a,old,older,t,oldt;
if(x<=1.)return -1.;
v= 3.*pow(.5*x,.666666);vi=1./v;
sum=0.;
term=1.;a=1.;old=0.;t=1.e10;
for(i=0;i<itmax;i++)
	{
	oldt=t;	
	t=term*a;
	if( abs(t)>abs(oldt)) break;
	sum+=t;
	term*=vi;
	older=old;
	old=a;
	if(!i) a=(3.*n*n+3.*n-1.)/12.;
	else
		{
		k=i+1;
		a=(.5*(n-2*k)*(2*k+3-n)*(2*k+3+2*n)*older
			-(12*k*k+36*k-3*n*n-3*n+25)*old)/(12*(k+2));
		}
	if(abs(t)<1.e-6|| abs(t)<1.e-5*abs(sum))break;
	itkt=i;
	}
return sqrt(pi/3.)*pow(v/3.,.5*n)*exp(-v)*sum;
}

double af1(x) double x;
{
if(x<3.)return af1s(x);
return afa(x,1.);
}
double af2(x) double x;
{
if(x<2.25)return af2s(x);
return afa(x,2.);
}
double af3(x) double x;
{
if(x<3.)return af3s(x);
return afa(x,3.);
}

double af(x,n) int n; double x;
{
int m;
double new,old,older,oldest;
if(n==1)return af1(x);
if(n==2)return af2(x);
if(n==3)return af3(x);
oldest=af3(x);older=af2(x);old=af1(x);

for(m=4;m<=n;m++)
	{
	new=.5*( (m-1)*older+x*oldest);	
	oldest=older;older=old;old=new;
	}
return new;

}

