/* functions defined in terms of Confluent Hypergeometric Function

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

Jbessel  I Associated Bessel function for complex order,argument
Ibessel  I Associated Bessel function for complex order,argument
Kbessel  K Associated Bessel function for complex order,argument
Airy	Ai CAVEAT: SEE TEXT!!!!!!!!!!!!!!!!!!!!!!
BiAiry	Bi
The following functions are included, but have not been benchmarked 
against any standard: 

bateman
Charlier
cunningham
Toronto

*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

/* multiply by i */
#define CTI(to,fm) {(to).x=-(fm).y;(to).y=(fm).x;}

extern int argrot;
/* some cpow() below could be ctreal()*/

Jbessel(order,arg, ans) struct complex *order,*arg,*ans;
	{struct complex a,b,c,d,iarg;int top=-1;
	CLET(a,*order);CTI(iarg,*arg);
	a.x += .5;/* a= order+.5*/
	CTREAL(b,a,2.);CLET(c,iarg);CTREAL(c,c,2.);
	c1f1(&a,&b,&c,top,ans);
	/*printc(ans);printf(" =M\n");*/
	CTREAL(b,iarg,-1.);
	cexp(&b,&a);
	CMULT(b,*ans,a);
	CTREAL(a,*arg,.5);
	cpow(&a,order,&c);CMULT(a,b,c);
	CLET(c,*order);c.x+=1.;
	cgamma(&c,&b,&d);
	CDIV((*ans),a,b);
	return 0;
	}

Ibessel(order,arg, ans) struct complex *order,*arg,*ans;
	{struct complex a,b,c,d;int top=-1;
	CLET(a,*order);
	a.x += .5;/* a= order+.5*/
	CTREAL(b,a,2.);CLET(c,(*arg));CTREAL(c,c,2.);
	c1f1(&a,&b,&c,top,ans);
	/*printc(ans);printf(" =M\n");*/
	CTREAL(b,*arg,-1.);
	cexp(&b,&a);
	CMULT(b,*ans,a);
	CTREAL(a,*arg,.5);
	cpow(&a,order,&c);CMULT(a,b,c);
	CLET(c,*order);c.x+=1.;
	cgamma(&c,&b,&d);
	CDIV((*ans),a,b);
	return 0;
	}

Kbessel(order,arg, ans) struct complex *order,*arg,*ans;
	{struct complex a,b,c;
	CLET(a,*order);
	a.x += .5;/* a= order+.5*/
	CTREAL(b,a,2.);CLET(c,(*arg));CTREAL(c,c,2.);
	cu(&a,&b,&c,ans);
	cexp(arg,&a);
	CDIV(b,(*ans),a);/*b=U exp(-x)*/
	CTREAL(a,*arg,2.);
	cpow(&a,order,&c);CMULT(a,b,c);/* a= Uexp()(2x)^n*/
	CTREAL((*ans),a,sqrt(pi));
	return 0;
	}

Airy(z,ans) struct complex *z,*ans;
	{struct complex a,b,c,d,z3h;double factor;double pow(),sqrt();
	int oldarg;
	oldarg=argrot;
	CMPLX(a,5./6.,0.);
	CMPLX(b,5./3.,0.);
	CMPLX(c,1.5,0.);
	factor= pow(3.,-5./6.)*pow(2.,2./3.)/sqrt(pi);
	
	cpow(z,&c,&d);CTREAL(z3h,d,4./3.);/*z3h= 4/3 z^(3/2)*/
	if(z->y==0.)
		{if( z->x <0.) {z3h.x=0.;argrot=1;}
		else z3h.y=0.;
		}
	/*printc(&z3h);printf(" U arg=4/3z^3/2 \n");*/
	cu(&a,&b,&z3h,&c);
	/*printc(&c);printf(" U, argrot=%d\n",argrot);*/
	CMULT(a,c,(*z));
	CTREAL(a,a,factor);
	CTREAL(z3h,z3h,-.5);
	cexp(&z3h,&b);
	argrot=oldarg;
	CMULT(*ans,b,a);return 0;
	}

BiAiry(z,ans) struct complex *z,*ans;
	{struct complex a,b,c,d,x,c1,c2,a1,a2,one;int oldrot;
	CMPLX(a,0., pi/6.);
	CMPLX(b,0., 2.*pi/3.);
	CMPLX(one,1.,0.);
	cexp(&a,&c1);
	CDIV(c2,one,c1);
	cexp(&b,&a1);
	CDIV(a2,one,a1);
	CMULT(x,a1,*z);	CLET(a1,x);
	CMULT(x,a2,*z);	CLET(a2,x);
	oldrot=argrot;argrot=0;
	Airy(&a1,&a);
	Airy(&a2,&b);
	CMULT(c,a,c1);CMULT(d,b,c2);CADD( *ans, c,d);
	argrot=oldrot;
	return 0;
	}


bateman(nu,x,ans) struct complex *nu,*x,*ans;
	{
	struct complex a,b,c,d;
	CLET(a,*nu);CTREAL(a,a,-.5);
	CMPLX(b,0.,0.);
	CTREAL(c,*x,2.);
	cu(&a,&b,&c,&d);
	cexp(x,&a);CDIV(b,d,a);
	CTREAL(a,*nu,.5);a.x+=1.;
	cgamma(&a,&c,&d);CDIV((*ans),b ,d);return 0;
	}

cunningham(n,m,x,ans)struct complex *n,*m,*x,*ans;
	{
	struct complex a,b,c,d,diff;
	CTREAL(b,*m,.5);CSUB(diff,b,*n);/* diff= m/2-n */
	CLET(c,*m); c.x+=1.;
	cu(&diff,&c,x,&d);
	CMPLX(a,1.,0.);CSUB(a,a,diff);/* a=1-diff= 1+n-m/2*/
	cgamma(&a,&b,&c);CDIV(c,d,b);
	CTI(b,diff);CTREAL(b,b,pi);CSUB(d,*x,b);CTREAL(d,d,-1.);
	cexp(&d,&b);CMULT((*ans),b,c);
	return 0;
	}

/* note that n is integer-Pochhammer symbol does not make sense otherwise*/
charlier(n,nu,x,ans) struct complex *nu,*x,*ans;int n;
	{
	struct complex a,b,c,d,e,f;int top=40;
	double gam;
	CMPLX(a,(double)(-n),0.);CADD(b,*nu,a);b.x+=1.;
	c1f1(&a,&b,x,-1,&d);
	cpochhammer(&b,n,&a);CMULT(f,a,d);
	gam=1./sqrt(gamma((double)(n+1)));
	CMPLX(b,-n*.5,0.);
	cpow(x,&b,&a); CTREAL(a,a,gam);
	CMULT((*ans),a,f);
	return 0;
	}

toronto(m,n,r,ans) struct complex *m,*n,*r,*ans;
	{struct complex mm,nn,rsq,a,b,c,d;int tol=50;
	CMULT(rsq,*r,*r);
	CLET(mm,*m);mm.x+=1.;CTREAL(mm,mm,.5);/* mm= .5(m+1) */
	CLET(nn,*n);nn.x+=1.;
	c1f1(&mm,&nn,&rsq,tol,&a);
	CTREAL(c,rsq,-1.);
	cexp(&c,&d);CMULT(c,a,d);/* c=Mexp(-r^2) */
	cgamma(&mm,&a,&b);CMULT(b,a,c);
	cgamma(&nn,&a,&c);CDIV(c,b,a);/*divide by n!=Gamma(n+1)*/
	CTREAL(a, *n,2.);CADD(a,a,*m); a.x-=1.;
	cpow( r,&a,&b); CMULT( (*ans),b,c);
	return 0;
	}

Laguerre( a, n, x,ans) struct complex *a,*x,*ans;int n;
	{
	struct complex p,q,r,s;
	CMPLX(q,n+1.,0.);
	cgamma(&q,&r,&p);/* r=Gamma(n+1)*/
	CMPLX(p,-n,0.);
	CLET(q,*a);q.x+=1.;/* q= a+1*/
	c1f1(&p,&q,x,-1,&s);
	CDIV( p,s,r);
	cpochhammer(&q,n,&s);
	CMULT((*ans),s,p);
	return 0;
	}