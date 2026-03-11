/*  

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

FILE *out;

main()
{double x;int i;
out=fopen("rich.out","w");
fprintf(out," 1\n");
for(i=1;i<100;i++)
	{
	x=i*.1;
	fprintf(out,"%le %le\n",x,ritchie(x));
	}
infinite_loop
	{
	printf(" enter x for Ritchie integral\n");
	scanf("%le",&x);printf(" ritchie(%le)=%le\n",x,ritchie(x));
	}
}