/* complex elliptic integrals of 1st,2nd kind

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

cef(m,z,e,f,sig) struct complex *e,*z,*f;double m,sig;
{
double m1,k,alpha, phi,psi,lambda,sinh(),sin(),cos(),cot,sinhp,cscp,
tan(),atan(),asin(),b1,b2,b3,mu,b,c,r1,r2,rp,sqrt(),cotl,b4,b5;
struct complex du1,du2;
m1=1.-m;
phi=z->x;psi=z->y;
cot= 1./tan(phi);cot*=cot;
cscp= 1./sin(phi);cscp*=cscp;sinhp=sinh(psi);sinhp*=sinhp;
b= -(cot+m*cscp*sinhp-m1);
c= -m1*cot;
lambda=b*b-4.*c ;
if(lambda<0.){fprintf(stderr," cef root would be complex\n");
	e->x=e->y=f->x=f->y=errorcode;return 1;}
r1=-b*.5;r2=.5*sqrt(lambda);
cotl=r1+r2; if(cotl<0.)cotl=r1-r2;
if(cotl<0.){fprintf(stderr," cef can't find positive root\n");
	e->x=e->y=f->x=f->y=errorcode;return 2;}
lambda= atan( 1./sqrt(cotl));
if(m==0.) {mu=pi*.5;if(abs(cotl-cot)<1.e-6)mu=0.;}
else { b4= (cotl/cot)-1.;
	if(b4<0.){
	fprintf(stderr,
	" warn mu imaginary %e %e %e %e\n",b4,cotl,cot,lambda);
	mu=0.;}
	else mu=atan(sqrt(abs(b4/m)));}
/*k=sqrt(m); alpha=asin(k);*/
b= sin(lambda);
c=sin(mu);r1=cos(mu);
b4=1.-m*b*b; b5= 1.-m1*b*b;
b1= m*b*cos(lambda)*c*c*sqrt(abs(b4));
b2=c*r1*b4*sqrt(abs(b5));
b3=r1*r1+m*b*b*c*c;
tef(lambda,m,sig,&b,&c);
/*printf("tef real: f e  %e %e\n",b,c);*/
tef(mu,m1,sig,&r1,&r2);
/*printf(" b1 b2 b3 %e %e %e\n",b1,b2,b3);*/
f->x =b; f->y=r1;
e->x=c;
e->y=r2-r1;
if(b4>=0.)e->x+=(b1/b3);
else e->y+=(b1/b3);
if(b5>=0.)e->x+=(b2/b3);
else e->y+=(b2/b3);
return 0;
}

czeta(m,z,ans,sig) double m,sig;struct complex *z,*ans;
{struct complex e,f,phi,term,ee,ff;double tanh(),am(),u;
if(m==0.){ans->x=0.;ans->y=0.;return 0;}
phi.x=pi*.5;phi.y=0.;
cef(m,&phi,&e,&f,sig);
cef(m,z,&ee,&ff,sig);
CDIV(phi,e,f);CMULT(term,phi,ff);
CSUB((*ans),ee,term);
return 0;
}

