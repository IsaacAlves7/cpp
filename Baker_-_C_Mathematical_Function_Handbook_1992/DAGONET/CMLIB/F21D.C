/* test driver 2F1

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#define cdigam cdigamma

main()
{struct complex a,b,c,d,e;float f,g,h,i,j,k,l,m,n,o,p;double z,digamma();
printf(" test driver for Gauss hypergeometric & related functions\n");
while(1)
	{printf(" enter real x for digamma\n");scanf("%e",&p);
	if(p==0.)break; z=p;  CMPLX(a,z,0.);cdigam(&a,&b);
	printf(" digamma=%le digam=%le\n",digamma(z),digam(z));
	printc(&b);printf("=cdigam\n");
	}

while(1)
{printf(" enter as complex numb a,b,c,x\n");
scanf("%e%e%e%e%e%e%e%e",&f,&g,&h,&i,&j,&k,&l,&m);
if(f==0. && g==0. && l==0. &&m==0. &&j==0. &&k==0. &&h==0. &&i==0.)break;
CMPLX(a,f,g);CMPLX(b,h,i);CMPLX(c,j,k); CMPLX(d,l,m);
cf21(&a,&b,&c,&d,&e);printf("answer= %e %e ",e.x,e.y);
}

while(1)
{printf(" Legendre P:enter as complex numb mu,nu,x\n");
scanf("%e%e%e%e%e%e",&f,&g,&h,&i,&j,&k);
if(f==0. && g==0. &&j==0. &&k==0. &&h==0. &&i==0.)break;
CMPLX(a,f,g);CMPLX(b,h,i);CMPLX(c,j,k);
cp(&c,&a,&b,&e);
printc(&e);printf("=answer\n");
}
while(1)
{printf(" Legendre Q:enter as complex numb mu,nu,x\n");
scanf("%e%e%e%e%e%e",&f,&g,&h,&i,&j,&k);
if(f==0. && g==0. &&j==0. &&k==0. &&h==0. &&i==0.)break;
CMPLX(a,f,g);CMPLX(b,h,i);CMPLX(c,j,k); 
cq(&c,&a,&b,&e);
printc(&e);printf("=answer\n");
}
}