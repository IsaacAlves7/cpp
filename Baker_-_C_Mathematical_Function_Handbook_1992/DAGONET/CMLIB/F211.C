
/* hypergeometric function complex arguments


from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define CTI(to,fm) {(to).x=-(fm).y;(to).y=(fm).x;}
#define eps 5.e-8
#define tol 1.e-7
#define DEBUG 0
#undef DEBUG

#define cdigam cdigamma

cpochhammer(x,n,ans) struct complex *x,*ans;int n;
{/* cannot use gamma(x+n)/gamma(x) as x may be 0 or neg. integer*/
struct complex product,factor,temp;int i;
if(!n){CMPLX(*ans,1.,0.);return 0;}
CLET(product,*x);
if(n==1){CSET(ans,product);return 0;}
for(i=1;i<n;i++)
	{
	CLET(factor,*x);factor.x+=(double)i;
	CMULT(temp,factor,product);
	CLET(product,temp);
	}
CSET(ans,product);return 0;
}

/*#define tol 1.e-6*/

cf21(a,b,c,x,ans)struct complex  *a,*b,*c,*x,*ans;
{/* complex hypergeometric function 2F1(a,b;c;x)*/
double gamma(),digamma();int m,i;struct complex s,z;
int aa,bb,cc,aint,bint,cint,top,l,mpl;
struct complex term,sum,ratio,factor,p,q,r,t,term1,term2,cl;
double count;
top=40;
#ifdef DEBUG
printf(" entered f21\n");
printc(a);printc(b);printc(c);printc(x);printf(" arg f21\n");
#endif


if(cabs(*x)<1.e-8 || cabs(*a)<1.e-8 || cabs(*b)<1.e-8)
	 {CMPLX(*ans,1.,0.);return 0 ;}
if(abs(a->y)<1.e-8)
	{aa=a->x;if((a->x-aa)==0.)aint=1;else aint=0;}else {aint=0;}
if(abs(b->y)<1.e-8)
	{bb=b->x;if((b->x-bb)==0.)bint=1;else bint=0;}else {bint=0;}
if(abs(c->y)<1.e-8)
	{cc=c->x;if((c->x-cc)==0.)cint=1;else cint=0;}else {cint=0;}
/*degenerate cases */
/* c,a or b negative integers- finite series n=0 to m, m=-a or b A&S 15.4.2*/
if(cint && cc<=0)
	{
	if (aint && aa<=0  && aa>cc)
		{
#ifdef DEBUG
		printf(" 2f1 degenerate series case\n");
#endif
		top=-aa;goto series2;}
	if (bint && bb<=0 && bb>cc )
		{
#ifdef DEBUG
		printf(" 2f1 degenerate series case\n");
#endif
		top=-bb;goto series2;}
	}
/* if c is 0 or neg. integer, returns F/Gamma(c) A&S 15.1.2*/
if(cint && cc<=0)
	{m=-cc;/*m>=0*/ count=1./ gamma(2.+m);/*count=m+1 !*/
	m++;
	cpochhammer(a,m,&p);cpochhammer(b,m,&q);CMULT(r,p,q);
	CTREAL(r,r,count);CMPLX(s,(double)m,0.);cpow(x,&s,&t);CMULT(s,r,t);
	CMPLX(term,(double)m,0.);CADD(term,term,*a);
	CMPLX(sum,(double)m,0.);CADD(sum,sum,*b);
	CMPLX(ratio,(double)m+1.,0.);
	cf21(&term,&sum,&ratio,x,&factor);
	CMULT((*ans),factor,s);
	return 2;
	} 
/* a or b negative integers, c not- polynomial finite series A&S 15.4.1*/
if((aint && aa<0 )||(bint && bb<0)) goto series;
/* c-a integer? */
CSUB(p,*c,*a);
if(abs(p.y)<1.e-7 && abs(p.x-((int)p.x))<1.e-7)
	{/*c-a integer*/
#ifdef DEBUG
printf(" c-a integer\n");
#endif
	if(p.x==0.){/*c-a=0,f=1/(1-z)^(-b)*/
		   CMPLX(q,1.,0.);CSUB(q,q,*x);CSUB(p,p,*b);/*p was 0 now -b*/
		   cpow(&q,&p,ans);return 0;}
	if(p.x<0.)/* c-a negative integer- get polynomial/(1-z)^power A&S 15.3.3*/
		{CMPLX(q,1.,0.);CSUB(q,q,*x);CSUB(r,p,*b);/*r= c-a-b*/
		cpow(&q,&p,&term1);CSUB(s,*c,*b);
		cf21(&p,&s,c,x,&term2);CMULT((*ans),term1,term2);
		return 0;}
	/*c-a positive integer: on to appropriate form*/
	}
#ifdef DEBUG
printc(c);printf(" c after int c-a \n");
#endif

/*x==1., Real(c-a-b)>0, c!= 0,-1,-2,...*/
#ifdef DEBUG
if(abs(x->y)<1.e-8 ) printf(" real argument\n");
if(abs(x->x-1.)<1.e-8)printf(" Re(z)==1\n");
if( (c->x)>(a->x+b->x))printf(" Re(c)>Re(a+b) \n");
#endif

if((abs(x->y)<1.e-8) && (abs(x->x-1.)<1.e-8)&&( (c->x)>(a->x+b->x)))
	{if(!cint || c->x>0.)/* do if not real, or not int, or positive*/
		{/*A&S 15.1.20*/
#ifdef DEBUG
printf(" 2f1 x==1\n");
#endif
		cgamma(c,&p,&ratio);
		CSUB(q,*c,*a);cgamma(&q,&r,&ratio);
		CDIV(ratio,p,r);/*ratio= gamma(c)/gamma(c-a)*/
		CSUB(q,*c,*b);cgamma(&q,&r,&s);
		CSUB(p,*c,*a);CSUB(p,p,*b);cgamma(&p,&s,&q);
		CDIV(p,s,r);CMULT((*ans),p,ratio);/* ans GcGc-a-b/Gc-a/Gc-b*/
		return 0;
		}
	}
/* Erdelyi(4) 15.3.9 as check on other stuff*/

/*
#ifdef DEBUG
printc(c);printf(" c before(4) \n");
#endif

CMPLX(r,1.,0.);CDIV(z,r,*x);CSUB(r,r,z);
if(cabs(r)<1. && cabs(*x)>1.)
	{printf(" Erdeyli (4) as check?\n");
	CMPLX(p,1.,0.);CADD(p,p,*a);CADD(p,p,*b); CSUB(p,p,*c);
	CMPLX(q,1.,0.);CADD(q,q,*c); CSUB(q,q,*a); CSUB(q,q,*b);
	if(!((q.x-(int)q.x)==0. && q.x<0.)&&!((p.x-(int)p.x)==0. && p.x<0.))
		{
#ifdef DEBUG
printf(" (4) OK\n");
#endif
		CMPLX(s,1.,0.);CADD(s,s,*a); CSUB(s,s,*c);
		cf21(a,&s,&p,&r,&t);
#ifdef DEBUG
printc(&t);printf(" first call returned\n");
#endif

		cpow(x,a,&s);CDIV(term1,t,s);
		cgamma(c,&r,&s);
		CSUB(t,*c,*a);CSUB(t,t,*b);cgamma(&t,&s,&p);
		CMULT(q,t,s);CMULT(t,term1,q);
		CSUB(p,*c,*a);cgamma(&p,&r,&s);
		CSUB(p,*c,*b);cgamma(&p,&q,&s);CMULT(p,r,q);CDIV(term1,t,p);
#ifdef DEBUG
printc(&term1);printf(" first term\n");
#endif
		CMPLX(q,1.,0.);CADD(q,q,*c); CSUB(q,q,*a); CSUB(q,q,*b);
		CMPLX(r,1.,0.);CDIV(z,r,*x);CSUB(r,r,z);
		CSUB(s,*c,*a);
		CMPLX(t,1.,0.);CSUB(t,t,*a);
		cf21(&s,&t,&q,&r,&p);
#ifdef DEBUG
printc(&p);printf(" 2nd call\n");
#endif
		CMPLX(q,1.,0.);CSUB(q,q,*x);CSUB(r,*c,*a);CSUB(r,r,*b);
		cpow(&q,&r,&s);CMULT(t,s,p);
		CSUB(q,*a,*c);cpow(x,&q,&s);CMULT(term2,s,t);
		cgamma(c,&r,&s);
		CADD(t,*b,*a);CSUB(t,t,*c);cgamma(&t,&s,&p);
		CMULT(q,t,s);CMULT(t,term2,q);
		cgamma(a,&r,&s);
		cgamma(b,&q,&s);CMULT(p,r,q);CDIV(term2,t,p);
#ifdef DEBUG
printc(&term2);printf(" 2nd term\n");
#endif
		CADD((*ans),term1,term2);
		return 0;
		}
	}
*/

CLET(q,*x);q.x-=1.;
if(cabs(q)<1. && x->x >.5){ f211(a,b,c,x,ans);return 0;}
else if(cabs(*x)<1.)/*series solution A&S 15.1.1*/
	{
	series:
	if(!(cint && cc<=0))/* c not 0 or negative integer*/
		{
		series2:
#ifdef DEBUG
printf(" 2f1 series\n");
#endif
		CLET(p,*a);CLET(q,*b);CLET(r,*c);
#ifdef DEBUG
printc(c);printc(&r);printf(" initial c term\n");
#endif
		for(i=0,term.x=sum.x=1.,term.y=sum.y=0.;i<top;i++)
			{
			count=(double)(i+1);
#ifdef DEBUG
printc(&p);printc(&q);printc(&r);
#endif
			CMULT(s,p,q);CDIV(ratio,s,r);
#ifdef DEBUG
printc(&ratio);printf(" p,q,r,ratio\n");
#endif
			ratio.x/=count;ratio.y/=count;
			CMULT(s,*x,ratio);
#ifdef DEBUG
printc(&s);
#endif
			CMULT(ratio,term,s);CLET(term,ratio);
#ifdef DEBUG
printc(&ratio);printf(" s, term\n");
#endif
			/*term *=(x*(a+i)*(b+i)/((c+i)*count));*/
			CADD(sum,sum,term);
			/*CDIV(ratio,term,sum);*/
#ifdef DEBUG
printc(&term);printc(&sum);printf(" term,sum 2f1 series\n");
#endif
			if( cabs(term)< tol*cabs(sum)){CSET(ans,sum);return 0;}
			p.x+=1.;q.x+=1.;r.x+=1.;
			}
		printf(" 2f1 tolerance not met series \n");
		CSET(ans,sum);return 0 ;
		}
	/* if c negative integer but a or b=n<m use A&S 15.1.2
	F is proportional to Gamma(c) which is infinite.*/
	CMPLX(*ans,errorcode,0.);
	return 0 ;/* c negative integer, no fixing*/
	}
else if(cabs(*x)>1.) {f21big(a,b,c,x,ans);return 0;}
if((cabs(*x)-1.)<1.e-7)/* on unit circle*/
	{if( a->x+b->x-c->x>= 1.){printf(" 2f1 divergent on unit circle\n");
		CMPLX(*ans,errorcode,0.);return 0 ;}
	/*divergent Erdelyi p.57 2.1.1 Vol.1*/
#ifdef DEBUG
	printf(" f21: |x|=1, conditionally convergent, attempting series\n");
#endif
	goto series;/*conditional convergence, try series*/
	}
/* never should get here? */
CMPLX(*ans,errorcode,errorcode);return 0;
}
