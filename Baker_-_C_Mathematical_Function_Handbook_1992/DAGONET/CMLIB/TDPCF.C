/*  test driver parabolic cylinder function and relatives

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

FILE *outp;

main()
{
double a,x,w,u,v,z;int i;
outp=fopen("pcf.out","w");fprintf(outp,"5\n");
printf(" enter a\n");scanf("%le",&a);
for(i=0; i<50;i++)
	{x=i/10.;
	u=upcf(a,x);v=vpcf(a,x);
	w=wpcf(a,x);
	z=wpcf(a,-x);
	fprintf(outp,"%le %le %le %le %le\n",x,u,v,w,z);
	}
infinite_loop
	{printf(" enter a, x. x>1000 to quit\n");
	scanf("%le%le",&a,&x);
	if(x>1000.)break;
	u=upcf(a,x);v=vpcf(a,x);
	printf(" U= %le, V %le\n",u,v);
	w=wpcf(a,x);
	z=wpcf(a,-x);
	printf(" W(%le,%le)=%le and for -x:%le\n",a,x,w,z);
	}
return 0;
}