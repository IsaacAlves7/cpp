/*  test driver U.c

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
main()
{float e,a,b,x;int i;
double exp(),aa,bb,xx,z,er,eps=1.e-6,u,up,uabx(),sqrt();
/* Ei test*/
for (i=1;i<=10;i++)
	{ xx= i;aa= exp(-xx)*uabx(1.,1.,xx,eps,&up);
/*uabx(1.,1.,-xx,eps,&up); xx<0 not allowed
	aa*=-exp(xx);                      */
	bb=aa*xx*exp(xx);
	printf(" x=%e E1=%e scaled=%e\n",xx,aa,bb);
	}
for (i=1;i<=10;i++)
	{ xx= i;aa= exp(-xx*xx)*uabx(.5,.5,xx*xx,eps,&up)/sqrt(pi);
	bb=aa*xx*exp(xx*xx);
	printf(" x=%e erfc=%e scaled=%e\n",xx,aa,bb);
	}
while(1)
{
printf("enter a,b,x,eps");scanf("%f%f%f%f",&a,&b,&x,&e);
aa=a;bb=b;xx=x;                          eps=e;
printf(" %e\n",uabx(aa,bb,xx,eps,&up));
}
}
