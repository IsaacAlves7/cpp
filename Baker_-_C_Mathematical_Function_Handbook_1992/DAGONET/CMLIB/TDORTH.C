/*  test driver orthogonal polynomials

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

main()
{
double a,b,x;float aa,bb,xx;int n;
while(1)
	{
	printf(" Jacobi enter a,b n x\n");
	scanf("%e%e%d%e",&aa,&bb,&n,&xx);
	a=aa;b=bb;x=xx;
	if(n<0||x<0.||a<0.||b<0.)break;
	printf(" %e\n",Pjacobi(a,b,n,x));
	}
while(1)
	{
	printf(" C gegenbauer enter a n x\n");
	scanf("%e%d%e",&aa,&n,&xx);
	a=aa;x=xx;
	if(n<0||x<0.||a<0.)break;
	printf(" %e\n",Cgegenbauer(a,n,x));
	}
while(1)
	{
	printf(" gen. Laguerre enter a n x\n");
	scanf("%e%d%e",&aa,&n,&xx);
	a=aa;x=xx;
	if(n<0||x<0.||a<0.)break;
	printf(" %e\n",laguerre(a,n,x));
	}
while(1)
	{
	printf(" P legendre enter n x\n");
	scanf("%d%e",&n,&xx);
	x=xx;
	if(n<0||x<0.)break;
	printf(" %e\n",Plegendre(n,x));
	}
while(1)
	{
	printf(" Tcheby enter n x\n");
	scanf("%d%e",&n,&xx);
	x=xx;
	if(n<0||x<0.)break;
	printf(" %e\n",Tcheby(n,x));
	}
while(1)
	{
	printf(" U Tcheby enter n x\n");
	scanf("%d%e",&n,&xx);
	x=xx;
	if(n<0||x<0.)break;
	printf(" %e\n",Ucheby(n,x));
	}
while(1)
	{
	printf(" Hermite enter n x\n");
	scanf("%d%e",&n,&xx);
	x=xx;
	if(n<0||x<0.)break;
	printf(" %e\n",Hermite(n,x));
	}

}