/*  test driver Bessel zero

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

main()
{
double order,z[20];
int i;
order=0.;
zerobes(order,10,z,1);
for(i=0;i<10;i++)printf(" %d zero for J %f = %e\n",i+1,order,z[i]);
zerobes(order,10,z,2);
for(i=0;i<10;i++)printf(" %d zero for Y %f = %e\n",i+1,order,z[i]);
zerobes(order,10,z,3);
for(i=0;i<10;i++)printf(" %d zero for J' %f = %e\n",i+1,order,z[i]);
zerobes(order,10,z,4);
for(i=0;i<10;i++)printf(" %d zero for Y' %f = %e\n",i+1,order,z[i]);
order=0.5;
zerobes(order,10,z,1);
for(i=0;i<10;i++)printf(" %d zero for J %f = %e\n",i+1,order,z[i]);
order=1.;
zerobes(order,10,z,1);
for(i=0;i<10;i++)printf(" %d zero for J %f = %e\n",i+1,order,z[i]);
return 0;
}
