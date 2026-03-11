/* hypergeometric function

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#define cdigam cdigamma
#define tol 1.e-7
#define DEBUG 1
#undef DEBUG
/* handles cases associated with series in (1-x) rather than x or 1/x*/
f211(a,b,c,x,ans) struct complex *a,*b,*c,*x,*ans;
	{double gamma(),digamma();
	int top=200,l,mpl,m,i;
	struct complex term,sum,ratio,factor,p,q,r,s,t,term1,term2,cl,z;
	double count;

	CADD(q,*a,*b);CSUB(p,*c,q);m=p.x;/*p=c-(a+b)*/
	if(abs(p.y)<1.e-7 && abs(p.x-m)<1.e-2)
		{/* Erdelyi 2.10(12)A&S 15.3.10-11 near pole */
#ifdef DEBUG
printf(" 2f1 Erdelyi(12) case A&S 10 or 11\n");
#endif
		if(m>=0)
		{
		CMPLX(z,1.,0.);CSUB(z,z,*x);/*z=1-x*/
		cgamma(c,&term1,&term2);/*term1=Gc*/
		CMPLX(term2,0.,0.);/* zero out for m==0 case*/
		if(m)
			{
			CLET(p,*a);CLET(q,*b);CMPLX(r,1.-m,0.);
			for(i=1,term.x=sum.x=1.,term.y=sum.y=0.;i<m;i++)
				{
				count=(double)(i);/*i=n*/
				CMULT(s,p,q);CDIV(ratio,s,r);
				ratio.x/=count;ratio.y/=count;
				CMULT(s,z,ratio);
				CMULT(ratio,term,s);CLET(term,ratio);
				CADD(sum,sum,term);
				if( cabs(term)< tol*cabs(sum))break;
				p.x+=1.;q.x+=1.;r.x+=1.;
				}
			CMPLX(p,gamma((double)m),0.);
			CMPLX(q,(double)m,0.);
			CADD(r,q,*a);CADD(s,q,*b);
			cgamma(&r,&t,&q);cgamma(&s,&q,&term2);
			CDIV(term2,p,t);/* Gm/Gm+a */
			CDIV(p,term2,q);/* p= Gm/Gm+a/Gm+b*/
			CMULT(term2,p,sum);/*term2= Sum*G/GG*/
			}
			/* infinite sum first find h[0]*/
			clog(&z,&cl);
			CMPLX(sum, digamma(1.+m)+digamma(1.),0.);
			CSUB(sum,sum,cl);
			CMPLX(q,(double)m,0.);CADD(q,q,*a);cdigam(&q,&p);
			CSUB(sum,sum,p);
			CMPLX(q,(double)m,0.);CADD(q,q,*b);cdigam(&q,&p);
			CSUB(sum,sum,p);/*sum=h[0]*/
			CLET(p,*a);CLET(q,*b);CMPLX(r,(double)(m+1),0.);
			p.x+=(double)m;
			q.x+=(double)m;
			CMPLX(term,1./gamma((double)(m+1)),0.);
			CMULT(sum,term,sum);/* sum now n=0 term in sum*/
			for(i=0;i<top;i++)
				{
				count=(double)(i+1);
				CMULT(s,p,q);CDIV(ratio,s,r);
				ratio.x/=count;ratio.y/=count;
				CMULT(s,z,ratio);
				CMULT(ratio,term,s);CLET(term,ratio);
				CMPLX(factor,
				  digamma(count+1.)+digamma(count+m+1),0.);
				CLET(s,*a);s.x+=(m+count);
				cdigam(&s,&t);
				CSUB(factor,factor,t);
				CLET(s,*b);s.x+=(m+count);
				cdigam(&s,&t);
				CSUB(factor,factor,t);
				CSUB(factor,factor,cl);/*factor=h*/
				CMULT(t,factor,term);
				CADD(sum,sum,t);
				if( cabs(t)< tol*cabs(sum))break;
				p.x+=1.;q.x+=1.;r.x+=1.;
				}
			CMPLX(q,(double)m,0.);
			cpow(&z,&q,&r);CMULT(p,sum,r);
			CTREAL(p,p,((m%2)?-1.:1));
			cgamma(a,&t,&q);CDIV(r,p,t);
			cgamma(b,&t,&q);CDIV(q,r,t);
			CADD(p,q,term2);CMULT((*ans),p,term1);
			return;
		}
		/* m negative integer  Erdelyi(14) A&S 12*/
		m=-m;
		CMPLX(z,1.,0.);CSUB(z,z,*x);/*z=1-x*/
		cgamma(c,&term1,&term2);/*term1=Gc*/
		CMPLX(term2,0.,0.);/* zero out for m==0 case*/
		if(m)
			{
			CLET(p,*a);p.x-=m;CLET(q,*b);q.x-=m;
			CMPLX(r,1.-m,0.);
			for(i=1,term.x=sum.x=1.,term.y=sum.y=0.;i<m;i++)
				{
				count=(double)(i);
				CMULT(s,p,q);CDIV(ratio,s,r);
				ratio.x/=count;ratio.y/=count;
				CMULT(s,z,ratio);
				CMULT(ratio,term,s);CLET(term,ratio);
				CADD(sum,sum,term);
				if( cabs(term)< tol*cabs(sum))break;
				p.x+=1.;q.x+=1.;r.x+=1.;
				}
			CMPLX(q,(double)(-m),0.);
			cpow(&z,&q,&p);
			CTREAL(p,p,gamma((double)m));
			cgamma(a,&t,&q);cgamma(b,&q,&term2);
			CDIV(term2,p,t);/* term2= (1-x)^-m Gm/Ga */
			CDIV(p,term2,q)
			CMULT(term2,p,sum);/*term2= finite term Sum*G/GG*/
			}
			/* infinite sum first find h[0]*/
			clog(&z,&cl);
			CMPLX(sum, digamma(1.+m)+digamma(1.),0.);
			CSUB(sum,sum,cl);
			CMPLX(q,(double)m,0.);CADD(q,q,*a);cdigam(&q,&p);
			CSUB(sum,sum,p);
			CMPLX(q,(double)m,0.);CADD(q,q,*b);cdigam(&q,&p);
			CADD(sum,sum,p);/*sum=h[0]*/
			CLET(p,*a);CLET(q,*b);CMPLX(r,(double)(m+1),0.);
			CMPLX(term,1./gamma((double)(m+1)),0.);
			CMULT(sum,term,sum);/* sum now n=0 term in sum*/
			for(i=0;i<top;i++)
				{
				count=(double)(i+1);
				CMULT(s,p,q);CDIV(ratio,s,r);
				ratio.x/=count;ratio.y/=count;
				CMULT(s,z,ratio);
				CMULT(ratio,term,s);CLET(term,ratio);
				CMPLX(factor,
				  digamma(count+1.)+digamma(count+m+1),0.);
				CLET(s,*a);s.x+=(count);
				cdigam(&s,&t);
				CSUB(factor,factor,t);
				CLET(s,*b);s.x+=(count);
				cdigam(&s,&t);
				CADD(factor,factor,t);
				CSUB(factor,factor,cl);/*factor=h*/
				CMULT(t,factor,term);
				CADD(sum,sum,t);
				if( cabs(t)< tol*cabs(sum))break;
				p.x+=1.;q.x+=1.;r.x+=1.;
				}
			CTREAL(sum,sum,((m%2)?-1.:1));
			CLET(p,*a);p.x-=m;
			cgamma(&p,&t,&q);CDIV(r,sum,t);
			CLET(p,*b);p.x-=m;
			cgamma(&p,&t,&q);CDIV(q,r,t);
			CADD(p,q,term2);CMULT((*ans),p,term1);
			return;
		}
	/* c not a+b +/- positive integer*/
	/* A&S 15.3.6  z=1.-x  A1 A2 coef. Erdelyi 2.10(1)*/
#ifdef DEBUG
printf(" 2f1 Erdelyi (1) A&S 6\n");
#endif
	cgamma(c,&term2,&ratio);
	CSUB(q,*c,*a);cgamma(&q,&r,&ratio);
	CDIV(ratio,term2,r);/*ratio= gamma(c)/gamma(c-a)*/
	CSUB(q,*c,*b);cgamma(&q,&r,&s);
	CSUB(p,*c,*a);CSUB(p,p,*b);cgamma(&p,&s,&q);
	CDIV(p,s,r);CMULT((*ans),p,ratio);/* ans= GcGc-a-b/Gc-a/Gc-b*/
	CMPLX(z,1.,0.);CSUB(z,z,*x);CMPLX(q,1.,0.);CADD(q,q,*a);CADD(q,q,*b);
	CSUB(q,q,*c);
	cf21(a,b,&q,&z,&t);CMULT(term1,*ans,t);/*first term*/

	cgamma(a,&r,&s);cgamma(b,&p,&s);CDIV(s,term2,r);CDIV(r,s,p);
	CADD(p,*a,*b);CSUB(p,p,*c);cgamma(&p,&s,&q);CMULT((*ans),s,r);
	/*ans= GcGa+b-c/GaGb*/
	CSUB(p,*c,*a);CSUB(q,*c,*b);CMPLX(r,1.,0.);CADD(r,r,*c);
	CSUB(r,r,*a); CSUB(r,r,*b);
	cf21(&p,&q,&r,&z,&term2);CMULT(p,term2,(*ans));/*p=A*f*/
	CSUB(q,*c,*a); CSUB(q,q,*b);
	cpow(&z,&q,&r);CMULT((*ans),r,p);CADD(*ans,*ans,term1);
	return;
	}
