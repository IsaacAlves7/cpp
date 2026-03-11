/*
Elliptic functions of Complex argument and modulus, using
Arithmetic-Geometric Mean (AGM)- Supplemental

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

Weiererstrass
InvWeier
ElExp
ElLog

*/
#include "Cmatrix.hpp"
#include <stdio.h>

#define maxterm 100
#define abstol 1.e-10
#define reltol 1.e-5
#define zerotol 1.e-14
#define errorcode -1.e60
#define pi 3.141592653589793238462643383279
#define lim 15.
#define max(a,b) ((a)>(b)?(a):(b))
#define infinite_loop for(;;)
#define rabs(x) ((x)<0.? -(x):(x))

complex tangent ( complex& x);
complex tanh(complex& x);
complex sinh(complex& x);
complex asin(complex&x);
complex acos(complex&x);
complex arctan(complex&x);
complex F( complex& phi0,complex& k);
void scdn( complex& ui,complex& m,complex& sn,complex& cn,complex& dn);

struct cmplx { double x,y;};
extern "C" ccubic(struct cmplx *,struct cmplx *,struct cmplx *,
	struct cmplx *,struct cmplx *,struct cmplx *);
extern "C" solvq(struct cmplx *,struct cmplx *,
		struct cmplx *,struct cmplx *);

void croots(complex& a,complex& b, complex& c,
	complex& b1, complex& b2, complex& b3)
	{
	struct cmplx a1,a2,a3,r1,r2, r3;
	a1.x=a.real();a1.y=a.imaginary();
	a2.x=b.real();a2.y=b.imaginary();
	a3.x=c.real();a3.y=c.imaginary();
	ccubic(&a1,&a2,&a3,&r1,&r2,&r3);
	b1=complex(r1.x,r1.y);
	b2=complex(r2.x,r2.y);
	b3=complex(r3.x,r3.y);
	}

void qroots(complex& a,complex& b,
	complex& b1, complex& b2)
	{
	struct cmplx a1,a2,r1,r2;
	a1.x=a.real();a1.y=a.imaginary();
	a2.x=b.real();a2.y=b.imaginary();
	solvq(&a1,&a2,&r1,&r2);
	b1=complex(r1.x,r1.y);
	b2=complex(r2.x,r2.y);
	}

//

complex Weierstrass(complex& x, complex& g2, complex& g3)
{
complex r1,r2,r3,a1,a2,a3,lambda,m,phi,k,beta,sn,cn,dn;
int rr=0;
a1=0.;a2=-.25*g2;a3=-.25*g3;
croots(a1,a2,a3,r1,r2,r3);
//r1.print("r1=","\n");r2.print("r2=","\n");r3.print("r3=","\n");
if( rabs(r2.imaginary()) < zerotol )
	{
	rr++;
	}
if( rabs(r1.imaginary()) < zerotol )
	{
	beta=r2; r2=r1;r1=beta;rr++;
	}
if( rabs(r3.imaginary()) < zerotol )
	{
	beta=r3; r3=r2;r2=beta;rr++;
	}
// below real m 1 real, 2 conj. roots r2=real root
//printf(" %d real roots\n",rr);
if(rr!= 1) goto alt;
beta=r2;
phi=( a2+(beta*(a1*2.+3.*beta) ))^.5;
lambda= phi^.5;
if(phi.abs() > zerotol )
	{
	m=(.5-.125*(6.*beta+2.*a1)/phi);
	if( rabs(m.imaginary()) < abstol ) m=complex( m.real(),0.);
	//m.print("m=","\n");
	x= -2.* x*lambda;
	//x.print("x=","\n");
	scdn(x,m,sn,cn,dn);
	// cn= cos phi. want x
	//cn.print(" cn=","\n");
	return (beta + phi*(1.+cn)/(1.-cn));
	}
alt:
//below real m  for 3 real roots or none
if( (x-r3).abs() > zerotol && (r1-r3).abs()>zerotol) goto ok;
else	
	{
	lambda=r1;r1=r3;r3=lambda;
	}
if( (x-r3).abs() > zerotol && (r1-r3).abs()>zerotol) goto ok;
else
	{
	lambda=r2;r2=r3;r3=lambda;
	}
if( (x-r3).abs() > zerotol && (r1-r3).abs()>zerotol) goto ok;
fprintf(stderr, " trouble in Weierstrass\n");
return complex(errorcode,errorcode);

ok:
// (for all real roots) r1>r2>r3
r1.print("r1=","\n");
r2.print("r2=","\n");
r3.print("r3=","\n");
if( r1.real()<r2.real())
	{lambda=r1;r1=r2;r2=lambda;}
if( r2.real()<r3.real())
	{lambda=r2;r2=r3;r3=lambda;}
if( r1.real()<r2.real())
	{lambda=r1;r1=r2;r2=lambda;}
//printf(" alternative method Weierstrass\n");
//r1.print("r1=","\n");
//r2.print("r2=","\n");
//r3.print("r3=","\n");
lambda= ((r1-r3)^.5) *.5;
//lambda.print(" lambda=","\n");
m= (r2-r3)/(r1-r3);
//m.print("alternative m=","\n");
k=m^.5;
phi= -x*lambda*2.;
scdn(phi,m,sn,cn,dn);
//k= beta + phi*(1.+cn)/(1.-cn);
lambda= r3+(r1-r3)/(sn*sn);//byrd and friedman 238
return lambda;
}

complex ElLog(complex& x, complex& a, complex& b)
{
complex r1,r2,r3,a1,a2,a3,lambda,m,phi,k;
a1=a;a2=b;a3=0.;
//croots(a1,a2,a3,r1,r2,r3);
r3=complex(0.,0.);
qroots(a,b,r1,r2);
if( (x-r3).abs() > zerotol && (r1-r3).abs()>zerotol) goto ok;
else
	{
	lambda=r1;r1=r3;r3=lambda;
	}

if( (x-r3).abs() > zerotol && (r1-r3).abs()>zerotol) goto ok;
else
	{
	lambda=r2;r2=r3;r3=lambda;
	}
if( (x-r3).abs() > zerotol && (r1-r3).abs()>zerotol) goto ok;
fprintf(stderr, " trouble in ElLog\n");
return complex(errorcode,errorcode);
ok:
phi= ((x-r1)/(x-r3))^.5;// cos phi
phi= acos(phi);
lambda= ((r1-r3)^.5) *.5;
m= (r2-r3)/(r1-r3); k= m^.5;
return (- F(phi,k)/lambda);
}


complex InvWeier(complex& x, complex& g2, complex& g3)
{
complex r1,r2,r3,a1,a2,a3,lambda,m,phi,beta,k,lou;

a1=0.;a2=-.25*g2;a3=-.25*g3;
croots(a1,a2,a3,r1,r2,r3);
//r1.print("r1=","\n");r2.print("r2=","\n");r3.print("r3=","\n");
if( (x-r3).abs() > 1.e-20) goto ok;
else
	{
	lambda=r1;r1=r3;r3=lambda;
	}

if( (x-r3).abs() > 1.e-20) goto ok;
else
	{
	lambda=r2;r2=r3;r3=lambda;
	}
if( (x-r3).abs() > 1.e-20) goto ok;
fprintf(stderr, " trouble in InvWeierstrass\n");
return complex(errorcode,errorcode);
ok:
//phi= ((x-r1)/(x-r3))^.5;// cos phi
//phi= acos(phi);
//lambda= ((r1-r3)^.5) *.5;
//m= (r2-r3)/(r1-r3);
beta=r2;
phi=( a2+(beta*(a1*2.+3.*beta) ))^.5;
//phi.print(" lambda^2=phi=","\n");
lambda= phi^.5;
m=(.5-.125*(6.*beta+2.*a1)/phi);
if( rabs(m.imaginary()) < abstol ) m=complex( m.real(),0.);
//m.print("m=","\n");
lou=(x-beta+phi);
//lou.print(" denom=","\n");
if( lou.abs()<zerotol)
	{
	return errorcode;
	}
phi= acos( (x-beta-phi)/lou);
//phi.print("phi=","\n");
k=m^.5;
return (.5* F(phi,k)/lambda);
}

complex ElExp(complex& Q, complex& a, complex& b)
{
complex r1,r2,r3,a1,a2,a3,lambda,m,x,sn,cn,dn;
a1=a;a2=b;a3=0.;
//croots(a1,a2,a3,r1,r2,r3);
r3=complex(0.,0.);
qroots(a,b,r1,r2);
if( (x-r3).abs() > 1.e-20) goto ok;
else
	{
	lambda=r1;r1=r3;r3=lambda;
	}

if( (x-r3).abs() > 1.e-20) goto ok;
else
	{
	lambda=r2;r2=r3;r3=lambda;
	}
if( (x-r3).abs() > 1.e-20) goto ok;
fprintf(stderr, " trouble in ElExp\n");
return complex(errorcode,errorcode);
ok:
//phi= ((r2-r1)/(x-r3))^.5;// cos phi
//phi= acos(phi);
lambda= ((r1-r3)^.5) *.5;
m= (r2-r3)/(r1-r3);
x= -Q*lambda;
scdn(x,m,sn,cn,dn);
//x= beta + phi*(1.+cn)/(1.-cn);
x= r3+(r1-r3)/(sn*sn);
return x;
}
