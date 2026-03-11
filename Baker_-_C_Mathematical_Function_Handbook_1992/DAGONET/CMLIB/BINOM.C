/* binomial coefficients
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

double binom(double n, double m)
{
double k;
if(m>n) return 0l;
/* ?  if m=3,n=2: 2!/3!-1! -1!=Gamma(0)=infinity=> 0 ans*/
if(m==0. || m==n)return 1l;
if(m==1.)return n;
if(n-m== 1.)return n;
k=n+1.;
return gamma(k)/(gamma(k-m)*gamma(m+1.));
/* slow for large n:*/
k=n-1.;
return binom(k,m)+binom(k,m-1.);
}
