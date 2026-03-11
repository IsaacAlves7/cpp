/* integrals of struve H0, H0/t, L0

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

double ModStruveI(x) double x;  integral of L from 0 to x
double StruveI(x) double x;     integral of H from 0 to x
double StruveIot(x) double x;   integral of H/t from 0 to x

*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define tol 1.e-8

double ModStruveI(x) double x;
{/* integral from 0 to x of L0(x)*/
double z,term,factor,sum;int i;
if(x==0.)return 0.;
z=x*x;sum=z*.5;term=z*z/36.;
for(i=5;i<10000;i+=2)
	{
	sum+=term;
	if( abs(term)<tol*abs(sum))goto fini;
	factor=i;
	term*= z*(factor-1.)/(factor*factor*(factor+1.));
	}
fprintf(stderr," ModStruveI could not meet error tolerance\n");
fini:return 2./pi*sum;
}

double StruveI(x) double x;
{/* integral from 0 to x of H0(x)*/
double z,term,factor,sum;int i;
if(x==0.)return 0.;
z=x*x;sum=z*.5;term=-z*z/36.;
for(i=5;i<10000;i+=2)
	{
	sum+=term;
	if( abs(term)<tol*abs(sum))goto fini;
	factor=i;
	term*= -z*(factor-1.)/(factor*factor*(factor+1.));
	}
fprintf(stderr," StruveI could not meet error tolerance\n");
fini:return 2./pi*sum;
}

double StruveIot(x) double x;
{/* integral from x to infinity of H0(x)/x*/
double z,term,factor,sum;int i;
if(x==0.)return pi/2.;
z=x*x;sum=x;term=-z*x/27.;
for(i=5;i<10000;i+=2)
	{
	sum+=term;
	if( abs(term)<tol*abs(sum))goto fini;
	factor=i;
	term*= -z*(factor-2.)/(factor*factor*(factor));
	}
fprintf(stderr," StruveI/t could not meet error tolerance\n");
fini:return pi/2.*(1.-4./(pi*pi)*sum);
}


