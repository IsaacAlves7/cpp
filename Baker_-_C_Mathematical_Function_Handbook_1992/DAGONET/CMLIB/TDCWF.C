/* test driver Coulomb Wave functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

FILE *outc;

main()
{
int l;float e,r;double eta,rho,f,g,fl,gl,fls,gls;

outc=fopen("cwf.plt","w");
fprintf(outc," 3\n");
printf(" enter eta\n");scanf("%le",&eta);
for(l=0;l<50;l++)
	{
	rho= .1*l;
	g=coulombg(eta,rho,0);
	f=coulombf(eta,rho,0);
	fprintf(outc,"%le %le %le\n",rho,f,g);
	}

while(1)
	{
	printf(" enter l>=0 eta rho ");
	scanf("%d%f%f",&l,&e,&r);eta=e;rho=r;
	if(l<0)break;
	printf(" result %e %e\n",coulombf(eta,rho,l),
	/*,coulombg(eta,rho,l)*/
	coulombg(eta,rho,l));
	cwfa(l,eta,rho,&fl,&gl,&fls,&gls);
	printf(" asmpt %e %e %e %e\n",fl,gl,fls,gls);
	}
}


