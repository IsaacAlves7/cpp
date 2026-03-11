
/* hypergeometric function complex arguments

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define cdigam cdigamma
#define tol 1.e-10
#define eps 5.e-8
#define DEBUG 1
#undef DEBUG

f21big(ain,bin,c,x,ans)struct complex *ain,*bin,*c,*x,*ans;
{double gamma(),digamma();
struct complex *a,*b;
int top=1000,i,l,m,mpl;
struct complex term,sum,ratio,factor,p,q,r,t,z,term1,term2,cl,s;
double count,rr;
#ifdef DEBUG
printf(" entered f21\n");
#endif
a=ain;b=bin;
CMPLX(p,1.,0.);CDIV(z,p ,(*x));/* z=1/x*/
CSUB(p,*b,*a);m=p.x;
if(abs(p.y)<1.e-7 && abs(p.x-m)<1.e-2)
	{/* near pole neg. integer 1-b+a or 1-a+b */
	if(m<0){m=-m;b=ain;a=bin;}
	CSUB(q,*c,*a); mpl=q.x;	l=mpl-m-1;
	if(abs(q.y)<1.e-7 && abs(q.x-mpl)<1.e-7 && q.x>0 && l>=0)
		{/* |z|>1, c-a,b-a integers b=a+m,c=a+m+l+1 m,l nonegative
		 use Erdelyi 2.10(9) pp,109-110 NOT in A&S */
		CMPLX(term2,0.,0.);/* zero out for m==0 case*/
#ifdef DEBUG
printf(" b-a integer, c-a also Erdelyi(9) m=%d l=%d\n",m,l);
#endif
		if(m)
			{
			CLET(p,*a);CMPLX(q,(double)(m),0.);CMPLX(r,m+l+1,0.);
			CMPLX(term,gamma((double)m)/gamma((double)(m+l+1)),0.);
			CLET(sum,term);
			CMPLX(t,0.,0.);CSUB(t,t,z);/*t=(-1/z)*/
			for(i=1;i<m;i++)
				{
				count=(double)(i);
				CMULT(s,p,q);CDIV(ratio,s,r);
				ratio.x/=count;ratio.y/=count;
				CMULT(s,z,ratio);
				CMULT(ratio,term,s);CLET(term,ratio);
				CADD(sum,sum,term);
				/*if( cabs(ratio)< tol*cabs(term))break;*/
				p.x+=1.;q.x+=1.;r.x+=1.;
				}
			cpow(&t,a,&r);
			CMULT(term2,r,sum);/*term2=1st finite sum*/
			}
#ifdef DEBUG
printc(&term2);printf(" m sum\n");
#endif
			/*infinite sum*/
			cpochhammer(a,m+l+1,&term);
			CTREAL(q,term,1./(gamma(2.+l+m)*gamma(2.+l)));
			CMPLX(p,(double)(l+1),0.);
			cpow(&z,&p,&t);/* t= z^-(l+1) */
			CMULT(term,q,z);
			CLET(sum,term);
			CMPLX(q,1.,0.);
			CMPLX(r,(double)(m+l+2),0.);
			CLET(p,*a);p.x+=(m+l+1);
			/* sum= n=l+1 term by here*/
			for(i=l+2;i<top;i++)
				{
				count=(double)(i);/*count=n*/
				CMULT(s,p,q);CDIV(ratio,s,r);
				ratio.x/=count;ratio.y/=count;
				CMULT(s,z,ratio);
				CMULT(ratio,term,s);CLET(term,ratio);
				CADD(sum,sum,term);
				if( cabs(term)< tol*cabs(sum))break;
				p.x+=1.;q.x+=1.;r.x+=1.;
				}
			CMPLX(s,0.,0.);CSUB(s,s,z);
			CLET(r,*a);r.x+=(double)m;
			cpow(&s,&r,&t);
			CMULT(term1,t,sum);
			CTREAL( term1,term1,((m+l+1)%2?-1.:1.));
#ifdef DEBUG
printc(&term1);printf(" infinite sum\n");
#endif
			CADD(term1,term1,term2);/*term1 first two terms*/
			/* last finite sum*/
			CMPLX(s,0.,0.);CSUB(s,s,*x);
			clog(&s,&cl);/*s= log(-x) */
			CMPLX(sum, digamma(1.+m)+digamma(1.)-digamma(1.+l),0.);
			CMPLX(q,(double)m,0.);CADD(q,q,*a);cdigam(&q,&p);
			CSUB(sum,sum,p);CADD(sum,sum,cl);/* sum=h'[0]*/
#ifdef DEBUG
printc(&sum);printf(" [log+] in n==0 term of lsum\n");
#endif
			cpochhammer(a,m,&t);
#ifdef DEBUG
printc(&t);printf(" pochhammer (a)m\n");
#endif
			CMPLX(s,(double)(-(m+l)),0.);
			cpochhammer(&s,m,&q);
#ifdef DEBUG
printc(&q);printf(" pochhammer -(m+l)m m=%d l=%d\n",m,l);
#endif
			CMULT(term,q,t);
			CTREAL(term,term,1./gamma((double)(m+1)));
#ifdef DEBUG
printf("  %e gamma in denom\n",gamma((double)(m+1)));
#endif
			CMULT(sum,term,sum);/*sum=(n=0 term)=term*[h'+log]*/
#ifdef DEBUG
printc(&sum);printf("  n==0 term of lsum\n");
#endif
			CLET(p,*a);p.x+=(double)m;
			CMPLX(q,(double)(-(l)) ,0.);
			CMPLX(r,(double)(m+1),0.);
			for(i=0;i<l;i++)/* n=i+1, sum to n=l */
				{
				count=(double)(i+1);
				CMULT(s,p,q);CDIV(ratio,s,r);
				ratio.x/=count;ratio.y/=count;
				CMULT(s,z,ratio);
				CMULT(ratio,term,s);CLET(term,ratio);
				CMPLX(factor,
				  digamma(count+1.)+digamma(count+m+1.)
				 -digamma(1.+l-count),0.);
				CLET(s,*b);s.x+=(count);/* b+n=a+m+n */
				cdigam(&s,&t);
				CSUB(factor,factor,t);
				CADD(factor,factor,cl);
				CMULT(t,factor,term);
				CADD(sum,sum,t);
				p.x+=1.;q.x+=1.;r.x+=1.;
				}
#ifdef DEBUG
printc(&sum);printf("  n==0 term of lsum\n");
#endif
			CMPLX(p,0.,0.);CSUB(p,p,z);
			cpow(&p,b,&r);CMULT(p,sum,r);
#ifdef DEBUG
printc(&p);printf(" scaled by (-x)^-b\n");
#endif
			CTREAL(p,p,1./gamma((double)(m+l+1)));
#ifdef DEBUG
printc(&p);printf(" l sum\n");
#endif
			CADD(term1,term1,p);/* all terms*/
			cgamma(c,&t,&r);CMULT(r,term1,t);/*c=a+m+l+1*/
			cgamma(b,&t,&q);CDIV((*ans),r,t);/*b=a+m*/
			return;
			/* end (9)*/
			}
	 /* b-a,integer c-a not 0 or pos. integer:
		Erdelyl 2.10(7) A&S 15.3.13(m=0) and 14*/
#ifdef DEBUG
printf(" Edelyi (7) b-a integer, c-a not 0 or pos integer,m=%d\n",m);
#endif
		CMPLX(term1,0.,0.);/* for m==0*/
		if(m)
			{
			CSUB(q,*c,*a);cgamma(&q,&p,&s);
			CMPLX(t,gamma((double)m),0.);
			CDIV(sum,t,p);/* sum=  n==0 term */
			CMPLX(term,1.,0.);CLET(q,*a);
			for(i=1;i<m;i++)/*i=n*/
				{
				count=(double)(i);/*n=count*/
				CMPLX(p,gamma((double)(m-i)),0.);
				CLET(r,*c);CSUB(r,r,*a);r.x-=count;
				cgamma(&r,&s,&t);/*s=Gc-a-n*/
				CDIV(factor,p,s);/*factor= Gm-n/Gc-a-n*/
				CTREAL(term,term,1./count);
				CMULT(s,z,term);CMULT(t,s,q);
				CLET(term,t);/*term=(a)n z^n/n!*/
				CMULT(ratio,factor,term);
				CADD(sum,sum,ratio);
				q.x+=1.;
				}
			CMPLX(p,0.,0.);CSUB(p,p,*x);cpow(&p,a,&r);
			CDIV(term1,sum,r);/*term1=first term*/
			}
#ifdef DEBUG
printc(&term1);printf("=finite term\n");
#endif
			/* infinite sum in 15.3.14*/
			CMPLX(s,0.,0.);CSUB(s,s,*x);
			clog(&s,&cl);
#ifdef DEBUG
printc(&z);printf("z=1/x\n");
#endif
			cpochhammer(a,m,&s);/* s= (a)m*/
#ifdef DEBUG
printc(&cl);printc(&s);printf(" log(-x)  poch (a)m\n");
#endif
			CMPLX(p,1.,0.);CADD(p,p,*a);CSUB(p,p,*c);/*p=1+a-c*/
			cpochhammer(&p,m,&t);
#ifdef DEBUG
printc(&t);printf(" poch 1+a-c)m\n");
#endif
			CMULT(term,t,s);CTREAL(term,term,1./gamma(1.+m));
#ifdef DEBUG
printc(&term);printf(" initial term value\n");
#endif
			CMPLX(factor,digamma(1.+m)+digamma(1.),0.);
			CADD(factor,factor,cl);
			cdigam(b,&s);CSUB(factor,factor,s);
			CSUB(s,*c,*b);
			cdigam(&s,&t);CSUB(factor,factor,t);
/*CMPLX(factor,1.,0.); old error-test to see effect*/
#ifdef DEBUG
printc(&factor);printf(" factor inn==0 term in infinite sum\n");
#endif
			CMULT(sum,term,factor);/*sum = n==0 term*/
			CMPLX(q,1.+m,0.);CSUB(q,q,*c);CADD(q,q,*a);
			/*CMPLX(r,(double)(m+1),0.);*/rr=(double)(m+1);
			CLET(p,*b);
#ifdef DEBUG
printc(&q);printc(&p);printf("=initial q,p\n");
#endif
			for(i=0;i<top;i++)
				{/*n=i+1=count*/
				count=(double)(i+1);
				CMULT(s,p,q);CTREAL(ratio,s,1./(rr*count));
				/*ratio.x/=count;ratio.y/=count;*/
				CMULT(s,z,ratio);
#ifdef DEBUG
printc(&s);printf(" ratio new term/old\n");
#endif
				CMULT(ratio,term,s);CLET(term,ratio);
				CMPLX(factor,
				  digamma(count+1.)+digamma(count+m+1),0.);
				CLET(s,*b);s.x+=(count);/*s=b+n=a+m+n*/
				if(abs(s.y)<eps){CMPLX(t,digamma(s.x),0.);}
				else cdigam(&s,&t);
				CSUB(factor,factor,t);
				CSUB(s,*c,*b);s.x-=count;/*s=c-b-n=c-a-m-n*/
				if(abs(s.y)<eps){CMPLX(t,digamma(s.x),0.);}
				else cdigam(&s,&t);
				CSUB(factor,factor,t);
				CADD(factor,factor,cl);/*factor=h+log*/
				CMULT(t,factor,term);
				CADD(sum,sum,t);
#ifdef DEBUG
printc(&sum);printc(&t);printc(&factor);printc(&term);printf(" sum,factor,term,i %d\n",i);
#endif
				if( cabs(t)< tol*cabs(sum))break;
				p.x+=1.;q.x+=1.;rr+=1.;
				}
			CTREAL(s,(*x),-1.);
			cpow(&s,b,&r);CDIV(p,sum,r);
			CSUB(q,*c,*a);
			cgamma(&q,&t,&s);CDIV(r,p,t);
			CADD(term1,term1,r);
#ifdef DEBUG
printc(&r);printc(&term1);printf(" inf sum tot\n");
#endif
			cgamma(c,&s,&t);CMULT(r,term1,s);
			cgamma(b,&s,&t);CDIV((*ans),r,s);
			return;
}
/* Erdelyi 2.10(2) A&S15.3.7  B1 & B2 */
#ifdef DEBUG
printf(" Erdelyi(2)\n");
#endif
cgamma(c,&term2,&ratio);
CSUB(q,*c,*a);cgamma(&q,&r,&ratio);
CDIV(ratio,term2,r);/*ratio= gamma(c)/gamma(c-a)*/
CSUB(q,*b,*a);cgamma(&q,&r,&s);
cgamma(b,&s,&q);
CDIV(p,r,s);CMULT(s,p,ratio);/* s=GcGb-a/Gc-a/Gb =B1*/
CMPLX(p,0.,0.);CSUB(p,p,*x);cpow(&p,a,&r);CDIV(t,s,r);/*t=(-x)^-aB1*/
CMPLX(p,1.,0.);CDIV(z,p,(*x));/* z=1/x*/
CMPLX(r,1.,0.);CADD(r,r,*a);CSUB(r,r,*b);/*r=1+a-b*/
CMPLX(q,1.,0.);CADD(q,q,*a);CSUB(q,q,*c);/*q=1+a-c*/
cf21(a,&q,&r,&z,&s);CMULT(term1,s,t);
/*term2= Gc at this point*/
cgamma(a,&r,&ratio);
CDIV(ratio,term2,r);/*ratio= Gc/Ga */
CSUB(q,*a,*b);cgamma(&q,&r,&s);
CSUB(p,*c,*b);cgamma(&p,&s,&q);
CDIV(p,r,s);CMULT(s,p,ratio);/* s=B2= GcGa-b/Ga/Gc-b*/
CMPLX(p,0.,0.);CSUB(p,p,*x);cpow(&p,b,&r);CDIV(t,s,r);/*t=B2*(-x)^-b*/
CMPLX(q,1.,0.);CADD(q,q,*b);CSUB(q,q,*c);
CMPLX(r,1.,0.);CADD(r,r,*b);CSUB(r,r,*a);
cf21(b,&q,&r,&z,&s);CMULT(term2,s,t);
CADD((*ans),term1,term2);
return;
}
