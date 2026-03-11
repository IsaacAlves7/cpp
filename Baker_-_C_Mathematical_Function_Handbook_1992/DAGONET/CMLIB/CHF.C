/* confluent hypergeometric function complex arguments

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

c1f1 confluent hypergeometric function of the first kind: M
cu      "            "           "      "  "  second kind: U

c1f1(a,c,x,top, ans).
	For top=-1, infinite series is used.
	in some applications, finite value of integer top needed
	If a is negative integer, polynomial.
	If c is negative integer, undefined unless
		a is a negative integer and |a|<|c|, i. e. a>c
	in this case, return 1F1/Gamma(c). function return -1

Cpow	x^n for complex x, n differs from cpow in that 
	while both use x^n= exp( n log x), for cpow
	log x has branch cut on neg. x axis, Cpow is adjustable.
	It adds two*pi*argrot, argot a global integer variable.	
	For argrot=1, argument of log x ranges from 0 to 2pi,
	rather than from -pi to pi.
	Cpow can therefore be used to find multiple roots by
	giving it argrot= 0,1,...
	For its use with cu, see text.

CAVEAT: YOU SHOULD SET argrot=0 BEFORE CALLING cu EXCEPT IN
VERY SPECIAL CIRCUMSTANCES (see text).
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"


/* to= i*fm*/
#define CTI(to,fm) {(to).x=-(fm).y;(to).y=(fm).x;}
#define cdigam cdigamma

#define tol 1.e-6
#define eps 1.e-7
#define topkt 1000

int argrot=0;/*default initialization. */

/* base ^ expon with ability to adjust log. principal value range*/

Cpow(struct complex *base,struct complex *expon,struct complex *ans)
{struct complex l,d;
clog(base,&l); l.y+=(2.*pi)*argrot; CMULT(d,l,*expon);cexp(&d,ans);
return 0;
}

c1f1(a,c,x,top,ans)struct complex  *a,*c,*x,*ans;int top;
{/* complex confluent hypergeometric function 1F1(a,b;c;x)*/
double oldterm,newterm,gamma(),digamma();int m,i;struct complex s,z;
int aa,bb,cc,aint,cint,l,mpl,topin;
struct complex term,sum,ratio,factor,p,q,r,t,term1,term2,cl;
double count;
struct complex savec,delta;
topin=top;
CLET(savec,*c);

if(top<0)top=topkt;
if(!top)
	{CMPLX((*ans),0.,0.);return 0;}
CLET(p,*a);CLET(r,*c);

if(cabs(*x)<eps){
	/*printf("small x=%le\n",cabs(*x));*/
	CMPLX((*ans),1.,0.);return 0;}
/* Kummer transformation- good idea?*/
/* suppress except for real arguments?*/
if(x->x <0.   /* &&  x->y ==0.*/ )
	{
	/*printf(" M:Kummer xfm\n");*/
	CTREAL(z,*x,-1.);  CSUB(p,*c,*a);
	c1f1(&p,c,&z,-1,&q);
	cexp(x,&p);
	CMULT((*ans),p,q);
	return 0;
	}

/* is a a negative integer or zero?*/
aint=0;
if( abs(a->y)<eps)
	{oldterm= a->x;aint=oldterm;
	if( aint <=0) aint= (abs(aint-oldterm))<eps ;
	else aint=0;
	}
/* is c a negative integer or zero?*/
cint=0;
if( abs(c->y)<eps)
	{oldterm= c->x;cint=oldterm;
	if( cint <=0) cint= (abs(cint-oldterm))<eps ;
	else cint=0;
	}
/*printf(" aint, cint top=%d %d %d\n",aint,cint,top);*/

if(cint && topin== -1 )
	{
	/*fprintf(stderr," warn: log case 1f1\n");*/
	if( !aint || cint>= aint)
		{
		CMPLX(s,1.,0.);CSUB(term,s,*c);CADD(sum,term,*a);
		c1f1(&sum,&term,x,top,&s);cpow(x,&term,&sum);
		CMULT(p,s,sum);CTREAL(p,p,1./gamma( 1.+term.x));
		cpochhammer(a,(int)(term.x),&s);
		CMULT(*ans,s,p);
		return -1;
		}
	}

if(cabs(*x)<=20.  || aint )
		{/*series*/
/*printf(" M:series\n");*/

		for(i=0,term.x=sum.x=1.,term.y=sum.y=0.;i<top;i++)
			{
			count=(double)(i+1);
			if(cabs(r)==0.)
				{
				fprintf(stderr," r=0 i=%d top=%d\n",i,top);
				/* once r== 0, will always be*/
				CSET(ans,sum); return 1;
				}
			CDIV(ratio,p,r);
			ratio.x/=count;ratio.y/=count;
			CMULT(s,*x,ratio);
			CMULT(ratio,term,s);CLET(term,ratio);
			CADD(sum,sum,term);
			/*CDIV(ratio,term,sum);*/
			if( cabs(term)< tol*cabs(sum))
				{
		/*printc(c);printf(" = c  series normal before cset\n");*/
				CSET(ans,sum);
				CSUB(delta,savec,*c);
				return 0;
				}
			p.x+=1.;r.x+=1.;
			}
if(top==topkt)fprintf(stderr," 1f1:requested accuracy not met\n");
/*otherwise, top was positive on entry for specified term count*/
/*printc(ans);printf("=M series\n");*/
		CSET(ans,sum);return 0;
		}
/*printf("M:asymptotic\n");*/
/* asymptotic expansion*/
CMPLX(p,1.,0.);CDIV(z,p,*x);
CSUB(p,*c,*a);CMPLX(r,1.,0.);CSUB(r,r,*a);
oldterm=1.e30;
for(i=0,term.x=sum.x=1.,term.y=sum.y=0.;i<top;i++)
		{
		count=(double)(i+1);
		CMULT(ratio,p,r);
		ratio.x/=count;ratio.y/=count;
		CMULT(s,z,ratio);
		CMULT(ratio,term,s);CLET(term,ratio);
		CADD(sum,sum,term);
		/*CDIV(ratio,term,sum);*/
		if( cabs(term)< tol*cabs(sum))break;
		newterm=cabs(term);
		if(oldterm>newterm)oldterm=newterm;
		else break;
		p.x+=1.;r.x+=1.;
		}
CSUB(q,*a,*c);cpow(x,&q,&p);CMULT(q,p,sum);
cexp(x,&p);CMULT(sum,p,q);
cgamma(a,&p,&q);CDIV(term1,sum,p);
CTREAL(z,z,-1.);
CLET(p,*a);CMPLX(r,1.,0.);CSUB(r,r,*c);CADD(r,r,*a);
oldterm=1.e30;
for(i=0,term.x=sum.x=1.,term.y=sum.y=0.;i<top;i++)
		{
		count=(double)(i+1);
		CMULT(ratio,p,r);
		ratio.x/=count;ratio.y/=count;
		CMULT(s,z,ratio);
		CMULT(ratio,term,s);CLET(term,ratio);
		/*term *=(x*(a+i)/((c+i)*count));*/
		CADD(sum,sum,term);
		/*CDIV(ratio,term,sum);*/
		newterm=cabs(term);
		if(oldterm>newterm)oldterm=newterm;
		else break;
		if( cabs(term)< tol*cabs(sum))break;
		p.x+=1.;r.x+=1.;
		}

cpow(x,a,&p);CDIV(q,sum,p);CLET(r,*a);CTREAL(r,r,-pi);CTI(s,r);
cexp(&s,&p);CMULT(sum,p,q);
CSUB(r,*c,*a);cgamma(&r,&p,&q);CDIV(term2,sum,p);
CADD(term2,term2,term1);
cgamma(c,&r,&s);CMULT((*ans),term2,r);
return 0;
}


cu(a,c,x,ans)struct complex  *a,*c,*x,*ans;
{double oldterm,newterm,digamma(),gamma();
struct complex z,s,term1,term2,term,sum; int m,i,top=100;
int aa,bb,cc,aint,l,mpl;
struct complex ratio,factor,p,q,r,t,cl;
double count;
if(cabs(*a)<eps)
	{CMPLX((*ans),1.,0.); return 0;}
if(cabs(*x)>8.)
	{/* asymptotic expansion*/
	/*printf(" asymptotic case\n");*/
	CMPLX(p,-1.,0.);CDIV(z,p,(*x));
	CLET(p,*a);CMPLX(r,1.,0.);CADD(r,r,*a);CSUB(r,r,*c);
	oldterm=1.e30;
	for(i=0,term.x=sum.x=1.,term.y=sum.y=0.;i<top;i++)
			{
			count=(double)(i+1);
			CMULT(ratio,p,r);
			ratio.x/=count;ratio.y/=count;
			CMULT(s,z,ratio);
			CMULT(ratio,term,s);CLET(term,ratio);
			/*term *=(x*(a+i)/((c+i)*count));*/
			CADD(sum,sum,term);
			/*CDIV(ratio,term,sum);*/
			if( cabs(term)< tol*cabs(sum))break;
			newterm=cabs(term);
			if(oldterm>newterm)oldterm=newterm;
			else break;
			p.x+=1.;r.x+=1.;
			}
		cpow(x,a,&p);CDIV((*ans),sum,p);return 0;
	}
m=c->x;
if( abs(c->y)<eps && abs(c->x-m)<eps)
	{/* c (b in A&S) is an integer*/
	/*printf(" integer case\n");*/
	
	if(m<1)
		{/*m =1-n  n=1-m>0  */
		/*printf(" integer case m<1 \n");*/
		CLET(p,*a);p.x+=1.-m; CMPLX(q,2.-m,0.);
		cu(&p,&q,x,&s);
		/*printc(&s);printf("= U(recursive)\n");*/
		CMPLX(p,1.-m,0.);cpow(x,&p,&r);CMULT((*ans),r,s);
		return 0;
		}
	/* m=n+1 m at least 1*/
	CMPLX(term1,1.,0.);
	CADD(term1,term1,*a);
	CSUB(term1,term1,*c);
	if( cabs(term1) < eps)
		{/* log. solution has pblm a=n for c=n+1. 
		use Kummer  13.1.29 */
		CMPLX(term2,2.,0.);CSUB(term2,term2,*c);
		cu(&term1,&term2,x,&p);
		CMPLX(term2,1.,0.);CSUB(term2,term2,*c);
		cpow(x,&term2,&term1);/* x^(1-b)*/		
		CMULT( *ans, term1,p);
		return 0;
		}	
	CMPLX(term1,0.,0.);
	if(m>1)
		{
		/*printf(" integer case m>1 \n");*/
		CMPLX(q,2.-m,0.);CLET(p,*a);p.x-=(m-1);
		c1f1(&p,&q, x,(m-2),&term1);/* was m-1*/
		/*printc(&term1);printf("= M m=%d\n",m);*/
		CTREAL(r,term1,gamma((double)(m-1)) );
		cgamma(a,&p,&q);CDIV(s,r,p);
		CMPLX(r,(m-1),0.);
		cpow(x,&r,&t);CDIV(term1,s,t);
		}
	clog(x,&s);
	c1f1(a,c,x,-1,&t);
	/*printc(&t);printf("= M\n");*/
	CMULT(term2,s,t);
	CLET(p,*a);CLET(r,*c);
	CMPLX(sum,-digamma(1.),0.);
	cdigam(a,&s);cdigam(c,&t);CADD(sum,sum,s);CSUB(sum,sum,t)
	for(i=0,term.x=1.,term.y=0.;i<top;i++)
			{
			count=(double)(i+1);
			CDIV(ratio,p,r);
			ratio.x/=count;ratio.y/=count;
			CMULT(s,*x,ratio);
			CMULT(ratio,term,s);CLET(term,ratio);
			CMPLX(factor,-digamma(count+m)-digamma(count+1.),0.);
			CLET(t,*a);t.x+=count;
			cdigam(&t ,&s);CADD(factor,factor,s);
			CMULT(t,ratio,factor);
			/*term *=(x*(a+i)/((c+i)*count));*/
			CADD(sum,sum,t);
			/*CDIV(ratio,term,sum);*/
			if( cabs(t)< tol*cabs(sum))break;
			p.x+=1.;r.x+=1.;
			}
	CADD(term2,term2,sum);
	CTREAL(term2,term2, ((m%2)?-1.:1.)/gamma((double)m) );
	CLET(p,*a);p.x-=(m-1);
	cgamma(&p,&q,&r);CDIV(s,term2,q);CADD((*ans),term1,s);
	return 0;
	}
/* x not large, b not integer*/
	{
	/*printf(" x not large, b no integer: use M\n");*/
	c1f1(a,c,x,-1,&s);
	cgamma(c,&t,&p);CDIV(q,s,t);
	CMPLX(p,1.,0.);CADD(p,p,*a);
	CSUB(p,p,*c);cgamma(&p,&s,&t);CDIV(term1,q,s);
	CMPLX(q,2.,0.);CSUB(q,q,*c);
	c1f1(&p,&q ,x,-1,&s);
	CMPLX(q,2.,0.);CSUB(q,q,*c);
	cgamma(a,&t,&p);CDIV(term2,s,t);
	cgamma(&q,&t,&s);CDIV(s,term2,t);
	CMPLX(p,1.,0.);CSUB(p,p,*c);
	/*CPOW(X,&p,&q);For Airy for Real arguments*/
	Cpow(x,&p,&q);
	CMULT(term2,s,q);
	CSUB(term1,term1,term2);CTREAL(term1,term1,pi);
	CTREAL(q,*c,pi);csin(&q,&s);
	CDIV((*ans),term1,s);
	return 0;
	}
}
