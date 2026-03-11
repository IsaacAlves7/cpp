/*  Toroidal P |x|>1

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

*/
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#define tol 1.e-8

double Ptoroidal(n,x) int n;double x;
{/*based on 7.10.8 of Lebedev */
double alpha,sum,fsum,term,factor,trm;int k,nmk;
if(x<0.)return errorcode;
alpha=arccosh(x);
fsum=0.;factor= exp(-2.*alpha);term=1.;
for(k=0;k<n;k++)
	{
	nmk=n-k;
	fsum+= gamma((double)nmk)*gamma(k+.5)
		/(gamma(nmk+.5)*gamma((double)(k+1)))*term;
	term*=factor;
	}
fsum*= exp( alpha*(n-.5))/pi;
trm=1.;
for(k=0,sum=0.;k<1000;k++)
	{
	nmk=n+k;
	term= gamma(nmk+.5)*gamma(k+.5)/(gamma(k+1.)*gamma(nmk+1.))
		*trm*(2.*alpha+digam(1.+k)-digam(.5+k)
		+digam(nmk+1.)-digam(nmk+.5));
	sum+=term;
	if(abs(term)<tol*abs(sum))goto fini;
	trm*=factor;
	}
fprintf(stderr," Ptoroidal did not meet tolerance\n");
fini:sum*= exp(-alpha*(.5+n))/(pi*pi);
return sum+fsum;
}

