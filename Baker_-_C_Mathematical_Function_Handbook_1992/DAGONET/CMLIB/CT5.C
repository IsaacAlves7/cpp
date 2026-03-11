/* Complex Elliptic Theta functions
	Copyright 1989,1991 Louis Baker. All rights reserved.
ctheta	theta functions complex v,q. v is (pi/2) times
	the value 2v used in Jahne & Emde's tables.
q	nome q given m for 0<= m,q <= 1
mq	m given q for 0<= m,q <= 1
emf	elliptic modular function m given q for complex q
emft	  "        "       "      m given t complex t, q= exp( i pi t)
amc	amplitude given complex m,u (see te.c for real values,
	set either by am() or in global double amplitude set by jef)

*/
#include "cmlib.h"
#include "protom.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define itmax 5000
#define tol 1.e-7

/* v is complex. q is complex, |q|<1, q= exp( -piK'/K)  v=pi*u/(2K)
where u is the argument of big-theta, sn, etc.*/

ctheta(v,q,ct1,ct2,ct3,ct4)struct complex *v,*q,*ct4,*ct3,*ct1,*ct2;
{/* t1,t2 converge more slowly that t3,t4 for larger |q|-
it may be desirable to split the evaluation into separate functions
in some applications if only t3 or t4 desired*/
struct complex sum,p2,d1,d2,one,term,logq,d3,d4,sum3;
struct complex sum1,sum2,qn1,s1,c1,ci,expiz;
double magv,si,old,older,oldest,worst;int vnonz;long int i;
old=older=oldest=10.; CMPLX(ci,0.,1.);
CMPLX(sum,1.,0.);CMPLX(d1,0.,2.);CMULT(p2,d1,(*v));
CMPLX(one,1.,0.);CLET(sum3,sum);
magv = cabs((*v));vnonz= magv>1.e-15;
if( vnonz)ctrig(v,&c1,&s1);
else {CMPLX(c1,1.,0.);CMPLX(s1,0.,0.);}
CLET(sum1,s1);CLET(sum2,c1);
CTREAL(sum1,sum1,2.);CTREAL(sum2,sum2,2.);CTREAL(c1,p2,.5);
if(vnonz)cexp(&c1,&expiz);
else{CLET(expiz,one);};
if(cabs((*q))>=1.)
	{
	ct1->x=ct2->x=ct3->x=ct4->x=errorcode;
	ct1->y=ct2->y=ct3->y=ct4->y=errorcode;return 1;
	}
si=-1.;clog(q,&logq);
for(i=1l;i<itmax;i++)
	{
	CTREAL(d1,logq,(double)i*i);cexp(&d1,&d2);
	CTREAL(d1,logq,(double)i*(i+1));cexp(&d1,&qn1);
	if(vnonz)
		{
		CTREAL(d1,p2,(double)i);cexp(&d1,&d3);
		if(cabs(d3)>1.e-20){CDIV(d4,one,d3);}
		else {CMPLX(d4,1.e20,0.);}
		CADD(d4,d4,d3);
		}
	else    {CMPLX(d3,1.,0.);CMPLX(d4,2.,0.);}
	CMULT(term,d2,d4);
/*printf(" cabs term %e %e\n",term.x,term.y);*/
	oldest=cabs(term);
	CADD(sum3,sum3,term);
	CTREAL(term,term,si);
	CADD(sum,sum,term);
	CMULT(d1,d3,expiz);
/*printf(" cabs d1 %e %e\n",d1.x,d1.y);*/
	if(cabs(d1)>1.e-20){CDIV(d4,one,d1);}
	else {CMPLX(d4,1.e20,0.);}
	CADD(c1,d1,d4);
	CMULT(term,qn1,c1);
	CADD(sum2,sum2,term);
/*printf(" cabs term %e %e\n",term.x,term.y);*/
	older=cabs(term);
	CSUB(term,d1,d4);CMULT(s1,term,ci);CTREAL(s1,s1,-1.);
	CMULT(term,qn1,s1);
	CTREAL(term,term,si);
	CADD(sum1,sum1,term);
/*printf(" cabs term %e %e\n",term.x,term.y);*/
	old=cabs(term);
/*	worst= max(old,max(older,oldest));*/
/*printf(" cabs term %e %e %e %e %e %e\n",
	sum1.x,sum1.y,sum2.x,sum2.y,sum.x,sum.y);*/
	if( cabs(sum1)*tol>old && cabs(sum2)*tol>older
	&& cabs(sum)*tol>oldest )break;
/*	oldest=older;
	older=old;*/
	si=-si;
	}
CSET(ct4,sum); CSET(ct3,sum3);
CTREAL(logq,logq,.25);cexp(&logq,&term);
CMULT( *ct1,term,sum1);CMULT(*ct2,term,sum2);
return 0;
}
/* given complex theta4, sn,cn,dn, theta1= sn*theta4(v)*theta2/theta3,
 2= cn*theta4(v)*theta2/theta4, 3= dn*theat4(v)*theta3/theta4*/

double q(m) double m;
{double k,alpha,cos(),c,asin(),src,eps,e4,q,a;
int flip;
/*k=sqrt(m); alpha= asin(k);     c=cos(alpha);*/
if(m==1.)return 1.;
/* better convergence for eps near 1/2 to switch k,k' and use
log q * log q' = pi^2*/
flip=0;
if(m> .5 )
	{flip=1;
	c=sqrt(m);
	}
else c= sqrt(1.-m);/* suffices for 0<=m<=1*/
src=sqrt(c);
eps=.5*(1.-src)/(1.+src);
e4= eps*eps; e4*=e4;
a=15.*e4;/* ratio of coeff. is abt 15 at final retained term,
asymptotes to 16. e4 is <= 1/16 . approx sum geometric series
for omitted terms as 1/(1-a)*/
q=eps*
 (1.+e4*(2.+e4*(15.+e4*(150.+e4*(1701.+e4*(20910.
  +e4*(268616.+e4*(3567400.+e4*(48555069.+e4*(673458874./(1.-a)))))))))));
if(flip)
	{q= exp(pi*pi/log(q)); }
return q;
}

#define itlimit 50
#define tol 3.e-8
int itermq;

double mq(q) double q;
{double eps,epsold,q8,q4,q2,nmult,num,denom,dmult,nterm,dterm,k;
q2=q*q;q4=q2*q2;q8=q4*q4;itermq=0;
nterm=q*q4*q4;num= q;nmult=q8*q8;
denom=1.;dterm=2.*q4;dmult=q4*q8;
epsold=0.;
while(itermq<itlimit)
	{num+=nterm;
	denom+=dterm;
	eps=num/denom;
	/*printf(" eps,epsold %le %le %le\n",eps,epsold,eps-epsold);*/
	if( abs(eps-epsold)<tol && eps<=.5)break;
	epsold=eps;
	nterm*=nmult;
	nmult*=q8;
	dterm*=dmult;
	dmult*=q8;
	itermq++;
	}
eps=min(.5,eps);
eps*=2.;
k = (1.-eps)/(1.+eps);
k*=k;/* epsold= (1+eps*(eps-2.))/(1+eps*(eps+2.));
printf(" k, alternate %le %le\n",k,epsold);*/
return 1.-(k*k);
}

emf(q,m) struct complex *q,*m;
{struct complex zero,rat,rat2,t1,t2,t3,t4;
CMPLX(zero,0.,0.);
ctheta(&zero,q,&t1,&t2,&t3,&t4);
CDIV(rat,t2,t3);CMULT(rat2,rat,rat);
CMULT((*m),rat2,rat2);
return 0;
}

emft(t,m)struct complex *t,*m;
{struct complex q,cip,arg;
CMPLX(cip,0.,pi);CMULT(arg,cip,(*t));cexp(&arg,&q);
emf(&q,m);
return 0;
}

amc(x,m,ans) double m;struct complex *x,*ans;
{
struct complex ci,z2,arg,sum,z,cn,dn;
cjef(x,m,&z,&cn,&dn);
CMULT(z2,z,z);
CMPLX(ci,0.,1.);CMPLX(arg,1.,0.);CSUB(arg,arg,z2);
clog(&arg,&sum);CTREAL(sum,sum,.5);cexp(&sum,&arg);
CMULT(sum,z,ci);CADD(sum,sum,arg);clog(&sum,&arg);
CMULT(sum,arg,ci);CTREAL((*ans),sum,-1.);
return 0;
}

/* for Am(z) use asin( sn(z)); if complex,
	asin = -i log(iz+sqrt(1-z^2)) */
