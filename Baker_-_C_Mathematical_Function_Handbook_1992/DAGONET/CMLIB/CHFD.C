/*

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

extern int argrot;

#define cdigam cdigamma
extern double digammin;

main()
{struct complex aa,bb,z,ans,w,f,g,h;double a,b,c,digamma(),d,e,ff;
int n,k;
argrot=0;
while(1)
	{printf(" enter x+iy for log\n");scanf("%le%le",&a,&b);
	if(a==0.&&b==0.)break;
	CMPLX(z,a,b);clog(&z,&ans);printc(&ans);printf(" log\n");
	}
while(1)
	{printf(" enter complex arg a,b z for confl. hypergom. 1f1 M\n");
	scanf("%le%le%le%le%le%le",&a,&b,&c,&d,&e,&ff);CMPLX(aa,a,b);
	if( a==0.&&b==0.)break;
	CMPLX(bb,c,d);CMPLX(z,e,ff);
	c1f1(&aa,&bb,&z,-1,&ans);
printc(&ans);printf("=answer\n");
	}
while(1)
	{printf(" enter complex arg a,b z for confl. hypergom. U\n");
	scanf("%le%le%le%le%le%le",&a,&b,&c,&d,&e,&ff);CMPLX(aa,a,b);
	if( a==0.&&b==0.)break;
	printf(" enter argrot\n");scanf("%d",&argrot);
	CMPLX(bb,c,d);CMPLX(z,e,ff);
	cu(&aa,&bb,&z,&ans);
printc(&ans);printf("=answer\n");
	}
	while(1)
		{
		printf(" enter digammin, x for digamma\n");scanf("%le%le",&digammin,&a);
		if(a==0.)break;
		printf(" psi(%le)=%le,floor=%le\n",a,digamma(a),digammin);
		}
	while(1)
		{
		printf(" enter digammin, y for digamma(1+iy)\n");scanf("%le%le",&digammin,&a);
		if(a==0.)break;
		CMPLX(w,1.,a);cdigam(&w,&f);
		printc(&f);
		}
	while(1)
		{
		printf(" K:enter order,arg\n");scanf("%le%le%le%le",&a,&b,&c,&d);
		if(a==0.&&b==0. && c==0. &&d==0.)break;
		CMPLX(w,a,b);
		CMPLX(g,c,d);
		Kbessel(&w,&g,&f);
		printc(&f);
		}
	while(1)
		{
		printf(" J:enter order,arg\n");scanf("%le%le%le%le",&a,&b,&c,&d);
		if(a==0.&&b==0. && c==0. &&d==0.)break;
		CMPLX(w,a,b);
		CMPLX(g,c,d);
		Jbessel(&w,&g,&f);
		printc(&f);
		}
	while(1)
		{
		printf(" I:enter order,arg\n");scanf("%le%le%le%le",&a,&b,&c,&d);
		if(a==0.&&b==0. && c==0. &&d==0.)break;
		CMPLX(w,a,b);
		CMPLX(g,c,d);
		Ibessel(&w,&g,&f);
		printc(&f);
		}
	while(1)
		{
		printf(" Airy:arg\n");
		/*scanf("%le%le%d",&a,&b,&argrot);*/
		scanf("%le%le",&a,&b);
		if(a==0.&&b==0.)break;argrot=0;
		/*printf(" argrot=%d \n",argrot);*/
		CMPLX(w,a,b);
		Airy(&w,&f);
		printc(&f);
		}
	while(1)
		{
		printf(" Bi(Airy):arg\n");
		scanf("%le%le",&a,&b);
		if(a==0.&&b==0.)break;
		CMPLX(w,a,b);argrot=0;
		BiAiry(&w,&f);
		printc(&f);
		}
	while(1)
		{
		printf(" gamma:arg\n");scanf("%le%le",&a,&b);
		if(a==0.&&b==0.)break;
		CMPLX(w,a,b);
		cgamma(&w,&f,&g);
		printc(&f);printc(&g);
		}
	while(1)
		{
		printf(" Erfc via Whittaker:arg\n");scanf("%le%le",&a,&b);
		if(a==0.&&b==0.)break;
		CMPLX(w,a,b);
		CMULT(h,w,w);CMPLX(aa,-.25,0.);CMPLX(bb,.25,0.);
		Wwhit(&aa,&bb,&h,&f);
		CTREAL(h,h,-.5);cexp(&h,&aa);CMULT(g,aa,f);
		/* Caveat: Whittaker & Watson p. 341
		definition of Erfc different from A&S 7.1.2
		*/
		CTREAL(g,g,1./sqrt(pi));
		ctreal(&w,-.5,&f);CMULT(aa ,f,g);
		printc(&aa);
		}
	while(1)
		{
		printf(" Poisson-Charlier: n,nu,x\n");scanf("%d%le%le",&n,&b,&a);
		/* A&C call x, nu what Erdelyi calls a, x respectively*/
		if(a==0.&&b==0. && !n)break;
		CMPLX(w,a,0.); CMPLX(f,b,0.);
		charlier(n,&f,&w,&g);
		printc(&g);printf("=from chf\n");
		ff=pow(a,(double)(-n*.5))*sqrt(gamma(1.+n));
		c=laguerre((b-n), n ,a);
		printf(" from Laguerre %le,coef %le ans: %le\n",c,ff,c*ff);
		/* check lag via recur*/
		c=1.;d=1.+b-n-a;  if(n==1)ff=d;if(!n)ff=c;
		for(k=1;k<n;k++)
			{
			ff= ((2*k+b-n+1.-a)*d-(k+b-n)*c)/(k+1.);
			c=d;d=ff;
			}
		printf(" Lag. by recursion=%le\n",ff);
		CMPLX(bb,b-n,0.);
		Laguerre(&bb,n,&w,&g);
		printc(&g);printf("=Lag. by 1F1\n");
		}
return 0;
}