/* Polylogarithm function

CAVEAT: DO NOT USE FOR DILOGARITHMS! WILL BE INACCURATE!
use Dlog or dilog instead.
(We do not automatically substitute a call for these when
n==1. to provide a crude cross-check on each)
Note also that n is not constrained to be an integer.

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define abs(x) ((x)<0.?-(x):(x))
extern int itkt;
#define tol 1.e-7

double polylog(x,n) double x,n;
{
int k,itmax=500;
double sum,term,newt,pow();
if(x< 0.) return errorcode;
sum=0.;
term=x;
for(k=1;k<itmax;k++)
	{
	newt= term/pow( (double)k,n);
	sum+=newt;
	/*printf(" k=%d newterm=%le sum=%le\n",k,newt,sum);*/
	if( abs(newt)< tol || abs(newt)< tol*sum)goto done;
	term*=x;
	}
itkt=itmax;
done:
itkt=k;
return sum;
}
