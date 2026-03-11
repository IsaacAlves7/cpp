/*  Legendre Q for complex arguments

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

/* to= i from */
#define CTI(to,from) {to.x=-from.y;to.y=from.x;}

cq(z,mu,nu,ans) struct complex *z,*mu,*nu,*ans;
{
double sqrt(),gamma(),x,y;
struct complex a,b,c,d,e,f,one,term1,term2,ccos,csin,p,q;
int j,k;
x=z->x;y=z->y;one.x=1.;one.y=0.;
if(abs(y)<1.e-8 && abs(x)<1.)
	{/* on cut  use Bateman 3.4(10) if mu not an integer*/
	/*printf(" on cut\n");*/
	if(!(abs(mu->y)<1.e-8 && abs(mu->x-((int)mu->x))<1.e-7))
	{
	/*printf(" munnot integer in effect A&S 8.1.6\n");*/
	CLET(d,one);CSUB(d,d,*z);CTREAL(d,d,.5);
	CMPLX(a,0.,0.);CSUB(a,a,*nu);
	CLET(b,one);CADD(b,b,*nu);
	CLET(c,one);CADD(c,c,*mu);
	cf21(&a,&b,&c,&d,&term1);
	CTREAL(d,*mu,.5);
	CMPLX(e, (1.-x)/(1.+x) ,0.);
	cpow(&e,&d,&f);CMULT(term2,term1,f);
	CTREAL(f,term2,.5);
	CLET(d,one);CADD(d,d,*nu); CADD(d,d,*mu);
	cgamma(&d,&e,&term1);CMULT(term1,e,f);
	CMPLX(d,0.,0.);CSUB(d,d,*mu);cgamma(&d,&e,&f);
	CMULT(f,term1,e);
	CLET(d,one); CADD(d,d,*nu);CSUB(d,d,*mu);
	cgamma(&d,&e,&term1);CDIV(term1,f,e);
/*printc(&term1);printf(" first term in sum\n");*/
	CLET(c,one);CSUB(c,c,*mu);
	cf21(&a,&b,&c,&d,&term2);
	CTREAL(d,*mu,.5);
	CMPLX(e, (1.+x)/(1.-x) ,0.);
	cpow(&e,&d,&f);CMULT(e,term2,f);
	CTREAL(e,e,.5);cgamma(mu,&a,&b);CMULT(f,a,e);
	CTREAL(a,*mu,pi);
	ctrig(&a,&b,&c);CMULT(term2,b,f);/* 2nd arg cos*/
	CADD( *ans,term1,term2);
/*printc(&term2);printf(" 2nd term in sum\n");*/
	return;
	}
/*printf(" integer mu on cut Erdelyi 3.4(12)\n");*/
	CLET(a,*nu);CADD(a,a,*mu);CTREAL(a,a,-.5);
	CLET(b,one);CADD(b,b,*nu);CSUB(b,b,*mu);CTREAL(b,b,.5);
	CLET(c,one);CTREAL(c,c,.5);/*c=1/2 b= (1+n-m)/2 a=-(n+m)/2*/
	CMULT(d,*z,*z);
	cf21(&a,&b,&c,&d,&term1);
/*printc(&term1);printf(" second f21\n");*/
	b.x+=.5;cgamma(&b,&e,&f);CDIV(f,term1,e);CTREAL(f,f,.5);
/*printc(&f);printf(" second term before tan/G((1-nu-mu)/2)\n");*/
	CADD(e,*mu,*nu);CTREAL(e,e,-.5);CTREAL(q,e,-pi);j=(int)e.x;
	ctrig(&q,&ccos,&csin);
	a.x+=.5;CLET(term2,f);
	CADD(p,*mu,*nu);k=(int)p.x;
	/* tan inf. if mu+nu odd integer*/
	if(!( abs(p.y)<1.e-8 && abs(p.x-(k))<1.e-8 && k%2))
		{
		cgamma(&a,&q,&term1);CDIV(term1,f,q);
		CMULT(f,term1,csin);CDIV(term2,f,ccos);
		}
	else
		{k>>=1;
		CTREAL(term2,term2,((k%2)?-1.:1.)*gamma(1.+k));
		}
/*printc(&term2);printf(" second term \n");*/

	CMPLX(c,1.5,0.);
	cf21(&a,&b,&c,&d,&term1);
	CMULT(f,*z,term1);
	CADD(a,a,*nu);cgamma(&a,&q,&c);CDIV(term1,f,q);
/*printc(&term1);printf(" first term before cot/G\n");*/
	if(!( abs(e.y)<1.e-8 && abs(e.x-(j))<1.e-8 && e.x<=1.e-8))
		{/*printc(&e);printf(" cut,NOT: mu+nu 0 or even integer\n");*/
		cgamma(&e,&p,&q);CDIV(q,f,p);CMULT(p,q,ccos);CDIV(term1,p,csin);
		}
	else
		{/*ratio of cot(pi*e/2)/Gamma(-e/2)*/
		/*printf(" mu+nu even integer %d\n",j);*/
		j=(-j);
		CTREAL(term1,term1,( (j)%2?-1.:1.)*gamma((double)j+1.));
		}
/*printc(&term1);printf(" first term \n");*/
	CSUB(term1,term1,term2);
	CTREAL((*ans),term1,sqrt(pi));
	if(cabs(*mu)>1.e-8)
		{
		CMPLX(c,2.,0.);cpow(&c,mu,&a);CMULT(term1,a,*ans);
		CSUB(d,one,d);CTREAL(a,*mu,-.5);cpow(&d,&a,&c);
		CMULT((*ans),term1,c);
		}
	return;
	}
/*printf(" off cut\n");*/
CLET(b,one);CADD(b,b,*mu); CADD(b,b,*nu); CTREAL(b,b,.5);
CMPLX(a,.5,0.);CADD(a,a,b);
CMPLX(f,1.,0.);CDIV(e,f,(*z));CMULT(f,e,e);/*f=1/z^2*/
CMPLX(c,1.5,0.);CADD(c,c,*nu);
/*printf(" cf call\n");printc(&a);printc(&b);printc(&c);printc(&f);*/
cf21(&a,&b,&c,&f,&e);/*printc(&e);printf(" answer\n");*/
CTREAL(e,e,sqrt(pi));
CMPLX(a,-1.,0.);CSUB(a,a,*mu);CSUB(a,a,*nu);cpow(z,&a,&b);CMULT(f,e,b);
CTREAL(b,*mu,.5);CMULT(a,*z,*z);CSUB(a,a,one);cpow(&a,&b,&c);CMULT(e,f,c);
CTI(a,(*mu));CTREAL(a,a,pi);cexp(&a,&b);CMULT(f,e,b);
CMPLX(a,-1.,0.);CSUB(a,a,*nu);CMPLX(b,2.,0.);cpow(&b,&a,&c);CMULT(e,f,c);
CLET(a,one);CADD(a,a,*mu);CADD(a,a,*nu);cgamma(&a,&c,&f);CMULT(f,e,c);
CMPLX(b,1.5,0.);CADD(b,b,*nu);cgamma(&b,&d,&e);CDIV((*ans),f,d);return;
}