/*  main() test driver for programs of Chapters 3-4
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <time.h>
#include <stdio.h>
#include "protom.h"
clock_t a,b,c;
long t1,t2;

#include "cmlib.h"
#include "complex.h"
#include "protom.h"

getcmplx( struct complex *z);


getcmplx(z) struct complex *z;
{double r,i;
scanf("%le%le",&r,&i); (z->x)=r;(z->y)=i;return 0;
}

extern double polyresid;

main()
{double x,y,c,d,e,f,r,s,t;int i;struct complex cz,cx,cy,cu,n,p;
/*a=clock();
for(i=0;i<10000;i++)
	y=tangnt(.1);
b=clock();
for(i=0;i<10000;i++)
	y=Tangnt(1.1);
c=clock();

printf(" using division(T): %ld  using first(t):%ld\n",c-b,b-a);
*/
printf(" test of cubic solver for imaginary coef\n");
CMPLX(p,0.,1.);CMPLX(n,0.,0.);ccubic(&n,&n,&p,&cx,&cy,&cz);
printc(&cx);printf("\n");
printc(&cy);printf("\n");printc(&cz);printf("\n");
printf(" test of cubic solver for imaginary coef\n");
CMPLX(p,0.,-1.);CMPLX(n,0.,0.);ccubic(&n,&n,&p,&cx,&cy,&cz);
printc(&cx);printf("\n");
printc(&cy);printf("\n");printc(&cz);
printf("\n test of quartic solver for imaginary coef\n");
CMPLX(p,0.,1.);CMPLX(n,0.,0.);cquartic(&n,&n,&n,&p,&cx,&cy,&cz,&cu);
printc(&cx);printc(&cy);printf("\n");
printc(&cz);printc(&cu);
printf("\n test of quartic solver for imaginary coef\n");
CMPLX(p,0.,-1.);CMPLX(n,0.,0.);cquartic(&n,&n,&n,&p,&cx,&cy,&cz,&cu);
printc(&cx);printc(&cy);printf("\n");
printc(&cz);printc(&cu);
printf("\n");
infinite_loop
	{printf(" enter real b,c,d cubic\n");
	scanf("%le%le%le",&c,&d,&e);
	if(c==0. && d==0. && e==0.) break;
	i=cubic(c,d,e,&r,&s,&t);
	printf(" ans components: %le %le %le %d resid %le\n",r,s,t,i,polyresid);
	}
infinite_loop
	{printf(" enter real b,c,d,e quartic\n");
	scanf("%le%le%le%le",&c,&d,&e,&f);
	if(c==0. && d==0. && e==0. && f==0.) break;
	i=quartic(c,d,e,f,&cx,&cy,&cz,&cu);
	printc(&cx);
	printc(&cy);
	printc(&cz);
	printc(&cu);
	printf(" quartic return value=%d resid=%le\n",i,polyresid);
	}
infinite_loop
	{printf(" power:enter x,y\n");scanf("%le%le",&x,&y);
	if(x>1000.)break;
	printf(" x^n=%le\n",power(x,y));
	i=y;
	printf(" x^(1/(int)n)=%le\n",root(x,i));
	}
infinite_loop
	{printf(" hypot:enter x,y\n");scanf("%le%le",&x,&y);
	if(x>1000.)break;
	printf(" Euclidean distance=%le\n",Euclidd(x,y));
	i=y;
	printf(" x^(1/(int)n)=%le\n",root(x,i));
	}
infinite_loop
	{printf(" cube rootenter x\n");scanf("%le",&x);
	if(x>1000.)break;
	printf("cbrt=%le\n",cube_rt(x));}
infinite_loop
	{printf(" tan enter x\n");scanf("%le",&x);if(x>10.)break;
	printf(" tan=%le %le \n",tangnt(x),tangent(x));}
infinite_loop
	{printf(" arc sin/cos enter x\n");scanf("%le",&x);if(x>1.)break;
	printf("asine=%le acos=%le\n",arc_sine(x),arc_cosine(x));}
infinite_loop
	{printf(" atan enter x\n");scanf("%le",&x);if(x==0.)break;
	printf("arc_tan=%le arc_tangent=%le\n",arc_tan(x),arc_tangent(2.*x,2.));}
infinite_loop
	{printf(" enter x\n");scanf("%le",&x);
	if(x==0.)break;
	printf("ln=%le\n",ln(x));}
infinite_loop
	{printf(" exp, hyper enter x\n");scanf("%le",&x);
	if(x>1000.)break;
	printf("exp=%le cosh %le sinh %le tanh %le\n",
		y=expon(x),hyper_cos(x),hyper_sin(x),hyper_tan(x));
	if(y==-errorcode || y==0.)break;
	}
infinite_loop
	{printf(" archyp enter x\n");scanf("%le",&x);
	if(x>1000.)break;
	printf("arc cosh =%le sinh %le tanh %le\n",
		arc_hyper_cos(x),arc_hyper_sin(x),arc_hyper_tan(x));
	if(y==-errorcode || y==0.)break;
	}
infinite_loop
	{printf(" sqrt enter x\n");scanf("%le",&x);
	if(x>1000.)break;
	if(x<0.)break;
	printf("sqrt=%le\n",square_rt(x));}
infinite_loop
	{printf("  enter complex z for clog\n");getcmplx(&cz);
	printc(&cz);
	if(cabs(cz)>100.)break;
	clog(&cz,&cx);
	printc(&cx);}

infinite_loop
	{printf("  enter complex z for ctrig\n");getcmplx(&cz);
	if(cabs(cz)>100.)break;
	ctrig(&cz,&cx,&cy);
	printc(&cx);}

infinite_loop
	{printf("  enter complex z for cexp\n");getcmplx(&cz);
	if(cabs(cz)>100.)break;
	cexp(&cz,&cx);
	printc(&cx);}
infinite_loop
	{printf("  enter complex z for ccosh\n");getcmplx(&cz);
	if(cabs(cz)>100.)break;
	ccosh(&cz,&cx);
	printc(&cx);}
infinite_loop
	{printf("  enter complex z for csinh\n");getcmplx(&cz);
	if(cabs(cz)>100.)break;
	csinh(&cz,&cx);
	printc(&cx);}
return 0;
}

