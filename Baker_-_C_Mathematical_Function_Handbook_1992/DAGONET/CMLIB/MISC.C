/*   sievert integral Ch27

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

extern int itkt,itmax;

double sievert(x,theta) double x,theta;
{
double sum,term,sq,en(),cos(),tt,y,z,kii();
int k,itmax=40;
z=cos(theta);
sq=z*z;
if(z!=0.) y=x/z;
else return errorcode;
sum=0.;
term=z;
for(k=0;k<itmax;k++)
	{
	tt=term*en(y,(double)((k+1)<<1) );
	sum+=tt;
	if(!k)term*=.5;
	else term*=((double)(2*k+1))/((double)((2*k+2)));
	term*=sq;
	if( abs(tt)<1.e-6 || abs(tt)<abs(sum)*5.e-7)break;
	}
z=kii(x);
/*printf(" sum=%e kii=%e\n",sum,z);*/
return z-sum;
}

