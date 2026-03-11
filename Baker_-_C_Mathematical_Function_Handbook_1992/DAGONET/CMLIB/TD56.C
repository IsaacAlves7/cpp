/*  test driver exponential integrals and gamma functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

extern struct complex ceiscaled;

main()
{float a,b,kk; int n;
double x,y,k,u,v;
struct complex c,d,e;
infinite_loop{
	printf(" E1: enter x<1000.\n");scanf("%le",&x);
	if(x>1000.)break;printf("e1=%le s=%le cf=%le\n",e1(x),e1s(x),en(x,1.));
	printf(" E1=%le\n",E1(x));
	}
infinite_loop{
	printf(" En: enter x<1000.,n\n");scanf("%le%d",&x,&n);
	if(x>1000.)break;
	printf(" En=%le\n",en(x,(double)n));
	}
infinite_loop{
	printf(" Ei: enter x<1000.\n");scanf("%le",&x);
	if(x>1000.)break;printf("Ei=%le\n",ei(x));
	}
infinite_loop{
	printf(" enter x<1000.,n  for alpha, beta\n");scanf("%le%d",&x,&n);
	if(x>1000.)break;printf("alpha=%le beta=%le\n",alpha(x,n),beta(x,n));
	}
infinite_loop{
	printf(" En complex arg:enter x,y k(==0. to end)");
	scanf("%e%e%e",&a,&b,&kk);
	if(kk==0.)break;x=a;y=b;k=kk;
	/* conventional En would have z^k-1 factor so we have additional z*/
	cei(x,y,k,1.e-5,&u,&v,&n);
	printf(" iter=%d ans= %e %e \n",n,u,v);
	CMPLX(c,x,y); cexpint(&c,k,1.e-5,&d,&n);
	printf(" En =");printc(&d);printf("\n");
	printf(" z^nExp(z)En=");printc(&ceiscaled);printf("\n");
	
	}
infinite_loop{
	printf(" gamma enter x<1000.\n");scanf("%le",&x);
	if(x>1000.)break;printf("gamma=%le %le\n",gamma(x),gammaqd(x));
	printf("digamma=%le %le\n",digamma(x),digam(x));
	printf("polygamma 1 and 2=%le %le\n",pg1(x),pg2(x));
	}
infinite_loop{
	printf(" complex gamma enter x,y<1000.\n");scanf("%le%le",&x,&y);
	if(x>1000.)break;CMPLX(c,x,y);cgamma(&c,&d,&e);
	printc(&d);printc(&e);
	printf(" gamma, log gamma \n");
	cdigamma(&c,&d);
	printc(&d);printc(&e);
	printf(" gamma, log gamma \n");
	}
infinite_loop{
	printf(" pochhammer enter x<1000,n\n");scanf("%le%d",&x,&n);
	if(x>1000.)break;printf("Pochhammer=%le\n",pochhammer(x,n));
	}
return 0;
}
