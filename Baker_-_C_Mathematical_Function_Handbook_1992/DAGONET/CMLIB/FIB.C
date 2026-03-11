/*  Fibonacci numbers
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
double fib(n) int n;
{
double a=1.618033988749894848204,b=-.618033988749894948204;
/*a= .5*(1+sqrt(5)), b=.5(1-sqrt(5))*/
double p,pow();
p=n;
return( pow(a,p)-pow(b,p))*.447213595;/*1/sqrt(5)*/
}
