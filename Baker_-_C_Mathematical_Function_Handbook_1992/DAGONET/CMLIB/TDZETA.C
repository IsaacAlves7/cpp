/*  test driver Riemann zeta(real), Bernoulli/Euler, 
		relatives.

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

main()
{
int i; double a,b,c,d;
infinite_loop
	{printf(" enter z for Riemann zeta =1 to end\n");
	scanf("%le",&a);
	if(a==1.)break;
	printf(" zeta=%le also %le\n",zeta(a), zeta2(a,1.));
	printf(" Catalan beta=%le lambda=%le eta=%le\n"
		,betacat(a),lambda(a),eta(a));
	}
infinite_loop
	{printf(" enter n for Bernoulli, Euler (-1 to end)\n");
	scanf("%d",&i);
	if(i==-1)break;
	printf(" Bernoulli= %le Euler=%le\n",bernoulli(i),euler(i));
	printf(" now enter x for poly\n");scanf("%le",&a);
	printf(" Bernoulli, Euler poly values %le %le\n"
		,bernpoly(i,a),eulerpoly(i,a));
	}
infinite_loop
	{printf(" enter z,mu for fermi integral -1 -1to end\n");
	scanf("%le%le",&a,&b);
	if(a==-1. && b==-1.)break;
	printf(" ifermi=%le\n",ifermi(b,a));
	}

}