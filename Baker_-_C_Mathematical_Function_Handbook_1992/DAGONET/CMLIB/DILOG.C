/*  Dilogarithm
Dlog differs from Spences's integral for n=2 given by dilog
in that Dlog(x)=dilog(1-x)

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define dl1 1.644934066848226
static double dcoef[17]=
{.005047192127203,.005300972587634,.004091615355944,.004815490327461,
.005966509196748,.006980881130380,.008260083434161,.009997129506220,
.012345919431569,.015625134938703,.020408155605916,.027777774308288,
.040000000124677,.062500000040762,.111111111110322,.249999999999859,1.};

double Dlog(x) double x;
{/* alternate version of Dilog- -Int from 0 to x of ln|1-y|/y*/
double f,u,y,z,log();int i;
if(x==1.)return dl1;
if(x==0.0)return 0.;
if(x>=2.)
	{z=1./x;f=-1.;u=log(x);u=-.5*u*u+3.289868133696453;}
else if(x>1.)
	{z=(x-1.)/x;f=1.;u=-.5*log(x)*log(z*x-z)+dl1;}
else if(x>.5)
	{z=1.-x;f=-1.;u=-log(z)*log(x)+dl1;}
else if(x>0.)
	{z=x;u=0.;f=1.;}
else if(x>=-1.)
	{z=x/(x-1.);f=-1.;u=log(1.-x);u=-.5*u*u;}
else
	{z=1./(1-x);f=1.;u=.5*log(z)*log(x*x*z)-dl1;}
y=.008048124718341*z+.008288074835108;
y=y*z-.001481786416153;
y=y*z-.000912777413024;
for(i=0;i<=16;i++)y=y*z+dcoef[i];
return f*y*z+u; 
}
