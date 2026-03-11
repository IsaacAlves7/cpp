/*  test driver for stirling numbers

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#include <stdio.h>

/* test driver stirling, fibonacci, binomial coef.

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

main()
{int n=1,k=1;
for(;;)
	{
	printf(" enter n,k ");
	scanf("%d%d",&n,&k);
	printf(" echo n, k %d %d\n",n,k);
	if(n<=0  || k<=0 )break;
	printf(" first %le %le second %le\n"
		,stirl1(n,k),stirlingf(n,k),stirl2(n,k));
	}
for(;;)
	{
	printf(" enter n for FIBONACCI ");
	scanf("%d",&n);if(n<=0)break;
	printf(" %le\n",fib(n));
	}
for(;;)
	{
	printf(" enter n,m for binomial coeff");
	scanf("%d%d",&n,&k);if(n<=0)break;
	printf(" %le\n",binom((double)n,(double)k));
	}
}