/*  Jacobian elliptic function complex argument 
	and Weierstrass P, sigma, zeta functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

cjef	complex jacobian elliptic functions sn,cn,dn
weier	Weierstrass function P,derivative P', eta=zeta(ometa),etc.
sigma	call weier and then use parameters to obtain Weierstrass sigma

ctld-logarithmic derivative of the theta functions, and
zeta, the Weirstrass zeta function, should be used with caution.  

Note zeta has a singularity at the origin

*/
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

cjef(u,m,sn,cn,dn) double m; struct complex *u,*sn,*cn,*dn;
{double w,s,c,d,s1,c1,d1;struct complex z;
jef(u->x,m,&s,&c,&d);
jef(u->y,1.-m,&s1,&c1,&d1);
w=1./( c1*c1+m*s*s*s1*s1);
CMPLX(z,w*(s*d1),w*(c*d*s1*c1));CSET(sn,z);
CMPLX(z,w*(c*c1),-w*(s*d*s1*d1));CSET(cn,z);
CMPLX(z,w*(d*c1*d1),-w*(s*m*s1*c));CSET(dn,z);
return 0;
}

static double ratiov;

double ratmp(m) double m;
{double x,y,z,p;
x=m-2.;y=2.*m-1.;z=m+1.;p=m*(m-1.)+1.;
if(p==0. || y==0.)
	{fprintf(stderr," zero denom. in ratmp m=%le x=%le y=%le\n",m,x,y);
	return errorcode;}
return .25*x*x*y*y*z*z/(p*p*p)-ratiov;
}

double ratmm(m) double m;
{double x,y,p;
/* return negative infinity for m=.5*/
if(m>=.5)return -1.e37;/* avoid overflw if single prec*/
x=32.*m*(m-1.)-1.;y=2.*m-1.;p=16.*m*(m-1.)+1.;
/*printf(" ratmm x y %le %le\n",x,y);*/
if(x==0. || y==0.)
	{fprintf(stderr," zero denom. in ratmm m=%le x=%le y=%le\n",m,x,y);
	return errorcode;}
return (p*p*p)/(x*x*y*y)-ratiov;
}

#define maxit 30
#define tol 1.e-5

double solvem(ratio,mu,ml,ratv) double ratio,mu,ml, (*ratv)();
{
 int i;
 double x,mbtm,mtop,value, topv,btmv,reltol=1.e-4,abstol=1.e-8,brent();
/* does it change sign from ml to mu? if not can't use brent*/
mbtm=ml;mtop=mu;     ratiov=ratio;
topv=ratv(mu);btmv=ratv(ml);
if(topv*btmv<0.)  return brent(mbtm,mtop,reltol,abstol,ratv);
/* cant use brent  use bisection????????*/
fprintf(stderr," solvem using bisection, not brent %le %le\n",topv,btmv);
 for(i=0;i<maxit;i++)
	{
	x=.5*(mbtm+mtop);
	value=ratv(x);
/*printf(" value=%le for x=%le\n",value,x);*/
	if(value==0.)return x;
	else if(value>0.) mbtm=x;
	else mtop=x;
	if( abs(value)<tol)break;
	}
value=.5*(mbtm+mtop);
/*printf(" solvem rv=%le %le %le %d\n",value, mbtm,mtop,i);*/
return value;

}

#define critm .066987398
/*critm=.5-sqrt(3)/4 */

double getm(g2,g3,d) int d; double g2,g3;
{double ml,mu,x,pow(),rat,rv;
if(d)
	{/* disc >0*/
	rat=g3*g3*27./(g2*g2*g2);
	ml=0.;mu=.5;
	return solvem(rat,mu,ml,ratmp);
	}
else
	{/* disc <0 */
	rat=(g2*g2*g2)/(27.*g3*g3);
	if (rat>=0.)
		{mu= critm ;
		ml=0.;
		}
	else
		{
		ml= critm ;
		mu=.5;
		}
/*printf(" getm: ml, mu %le %le\n",ml,mu);*/
	rv=solvem(rat,mu,ml,ratmm);
/*printf(" solvem returned %le\n",rv);*/
	return rv;
	}
return 0.;
}

weier(z,g2,g3,p,pp,mp,kp,ep,omegap,eta,ehp)
double g2,g3,*mp,*kp,*ep,*omegap,*ehp,*eta; struct complex *z,*p,*pp;
{
int dpos;
struct complex ci,x,q,qp,aux,sn,cn,dn,aux2;
double pow(),discr,d1,d2,sqrt(),k,e,m,omega,om2,e1,e2,e3,h2,c,ee,getm();

CMPLX(ci,0.,1.);
if(g3<0.) {  CMULT(x,ci,*z);weier(&x,g2,g3,&q,&qp,mp,kp,ep,omegap,eta,ehp);
	CMPLX(aux,0.,0.);CSUB(*p,aux,q);CSUB(*pp,aux,qp);
	return 1;
	}
d1= 27.*g3*g3;d2=(g2*g2*g2);
discr=d2-d1;
/*printf(" discr=%le\n",discr);*/
if( discr==0.) /* DEGENERATE CASE not doubly-periodic */
	{ /* g3>=0*/
	if( g3==0. /* && g2==0. */)
		{
		CMPLX(aux,1.,0.);CDIV(q,aux,*z);CDIV((*p),aux,*z);
		CDIV((*pp),*p,*z);CTREAL(*pp,*pp,-2.);return 0;
		}
	c=sqrt(abs(g2/12.));
	CMPLX(aux,c,0.); if(g2<0.){CMULT(q,ci,aux);CLET(aux,q);}
	/*p= -c+3c/ {sin [sqrt(3c)z]}^2*/
	ee=sqrt(3.*c);
	CTREAL( aux, *z,ee);csin(&aux,&aux2);
	CMULT(aux,aux2,aux2); CMPLX(aux2,3.*c,0.);CDIV(q,aux2,aux);
	CLET(*p,q);p->x -=c;
	ccot(&aux,&qp);CTREAL(qp,qp,-2.*ee);CMULT(*pp,qp,q);
	return 1;
	}
else if(discr>0.)
	{
	/*find 0<m<.5*/
	if(g3==0.)m=.5;
	else m=getm(g2,g3,1);
	tek(0,m,&k,&ee);
	omega= sqrt(sqrt(4./3.*(m*(m-1.)+1.)/g2))*k;
	e= k*k/(3.*omega*omega);
	e1=e*(2.-m);e3=-e*(1.+m);e=sqrt(abs(e1-e3));
/*printf(" omega=%e,e1=%e,e3=%e, sqrt(dif) %e\n",omega,e1,e3,e);*/
	*omegap=omega; *ehp=e;
	*eta= k/(3.*omega)*(3.*ee+(m-2.)*k);
	if((e1-e3)>0.){CMPLX(aux,e,0.);}
	else	{CMPLX(aux,0.,e);}
	CMULT(x,*z,aux);
	cjef(&x,m,&sn,&cn,&dn);
	CMULT(aux2,sn,sn);
	CMPLX(q,(e1-e3),0.);
	CDIV(qp,q,aux2);
	qp.x+=e3;CSET(p,qp);
	CMULT(q,aux2,sn);
	CDIV(qp,cn,q);
	CMULT(q,qp,dn);
	CTREAL(q,q,-2.*(e1-e3));
	CMULT(*pp ,q,aux);
	}
else    /* discr<0.*/
	{
	/* find m*/
	/*	printf(" discr<0\n");*/
	if(g2==0.) m= critm ;
	else m=getm(g2,g3,0);
/*printf(" after getm m=%le g2,3 %le %le\n",m,g2,g3);*/
	tek(0,m,&k,&ee);
	om2= k*pow( 8.*(2.*m-1.)*(32.*m*(m-1.)-1.)/(27.*g3),1./6.) ;
	if(om2==0.)
		{fprintf(stderr," om2==0. in weier\n");
		p->x=errorcode;p->y=errorcode;
		pp->x=errorcode;pp->y=errorcode;
		return 1;
		}
	e=k*k/(3.*om2*om2);
	e2=2.*(1.-2.*m)*e;
	h2= .75*e2/(.5-m);
	e=sqrt(abs(h2));
/*printf(" omega2=%e,e2=%e,h2=%e, sqrt(dif) %e\n",om2,e2,h2,e);*/
	*omegap=om2;*ehp=h2;
	*eta=k/(3.*om2)*(6.*ee+(4.*m-5.)*k);
	if(h2>0.){CMPLX(aux,e,0.);}
	else{CMPLX(aux,0.,e)};
	CMULT(x,*z,aux);
	CTREAL(x,x,2.);
	cjef(&x,m,&sn,&cn,&dn);
	CLET(q,cn);q.x+=1.;
	CMPLX(aux2,1.,0.);CSUB(aux2,aux2,cn);
	if(cabs(aux2)==0.)
		{fprintf(stderr," cn= %e %e divide error\n",cn.x,cn.y);
		p->x=errorcode;p->y=errorcode;
		pp->x=errorcode;pp->y=errorcode;
		return 1;
		}
	CDIV(qp,q,aux2);
	CTREAL(qp,qp,h2);qp.x+=e2;CSET(p,qp);
	CMULT(qp,aux2,aux2);
	if(cabs(qp)==0.)
		{fprintf(stderr," qp==0. in weier\n");
		p->x=errorcode;p->y=errorcode;
		pp->x=errorcode;pp->y=errorcode;
		return 1;
		}
	CMULT(q,sn,dn);CDIV(aux2,q,qp);
	CMULT(q,aux2,aux);
	CTREAL(q,q,-4.*h2);CSET(pp,q);
	}
*mp=m;*ep=e;*kp=k;
return 0;
}

ctld(v,q,t1,t2,t3,t4) struct complex *v,*q,*t1,*t2,*t3,*t4;
{/* log derv. theta functions*/
double si;
struct complex sin,cos,denom,qn,q2n,one,ti,fact1,fact2,q2,x,e2ui,e2uii;
struct complex q2e2,q2e2i,qe2,qe2i;
struct complex term1,term2,term3,term4,sum1,sum2,sum3,sum4;
int n,itmax=1000;
CMPLX(one,1.,0.); CMPLX(ti,0.,2.);
ctrig(v,&cos,&sin);
CDIV((*t1),cos,sin);CDIV((*t2),sin,cos);CTREAL(*t2,*t2,-1.);
CMPLX(*t4,0.,0.); CSET(t3,*t4);
/*printf(" cot %e %e\n",t1->x,t1->y);*/
CLET(qn,*q);CMULT(q2,qn,qn);CLET(q2n,q2);
CMPLX(sum1,0.,0.); CLET(sum2,sum1);CLET(sum3,sum1);CLET(sum4,sum1);si=-1.;
CMULT(x,ti,*v);cexp(&x,&e2ui);CDIV(e2uii,one,e2ui);
CSUB(x,e2uii,e2ui); /*if( cabs(x)<1.e-10) return;*/
CMULT(q2e2,q2n,e2ui);CMULT(q2e2i,q2n,e2uii);
CMULT(qe2,qn,e2ui);CMULT(qe2i,qn,e2uii);
if(cabs(qn)>=1.)printf(" warn ctld q= %e %e |q|>1\n",qn.x,qn.y);
for(n=1;n<itmax;n++)
	{
	CSUB(denom,one,q2n);
	CSUB(x,q2e2,q2e2i);CDIV(term1,x,denom);
	CLET(term2,term1);
	CTREAL(term2,term2,si);
	CSUB(x,qe2,qe2i);CDIV(term4,x,denom);
	CLET(term3,term4);
	CTREAL(term3,term3,si);
	CADD(sum1,sum1,term1);
/*printf(" sum1, term1 %e %e %e %e %d\n",sum1.x,sum1.y,term1.x,term1.y,n);*/
	CADD(sum2,sum2,term2);
	CADD(sum3,sum3,term3);
	CADD(sum3,sum4,term4);
	if(cabs(sin)>1.e-3){
	if( cabs(term4)<tol*cabs(sum2) && cabs(term1)<tol*cabs(sum1))break;}
	if( cabs(sum1)<tol && n>10)break;/* for u=0,pi/2 type cases*/
	CMULT(sin,q2e2,q2);CMULT(q2e2,sin,e2ui);
	CMULT(sin,q2e2i,q2);CMULT(q2e2i,sin,e2uii);
	CMULT(sin,qe2,qn);CMULT(qe2,sin,e2ui);
	CMULT(sin,qe2i,qn);CMULT(qe2i,sin,e2uii);
	si=-si;
	/* avoid an overflow*/
	if( cabs(term4)>1.e10 || cabs(term1)>1.e10)break;
	}
CMULT(sin,sum1,ti);
/*printf(" sum1 %e %e,term1 %e %e %d\n",sum1.x,sum1.y,term1.x,term1.y,n);*/
CMULT(cos,sum2,ti);
CMULT(*t3,sum3,ti);
CMULT(*t4,sum4,ti);
CADD(*t1,*t1,sin);CADD(*t2,*t2,cos);
return 0;
}

sigma(z,ans,k,omega,eta,m,g2,g3) struct complex *z,*ans;
double k,omega,eta,m,g2,g3;
{struct complex v,dum1,dum2,q,t1,t2,t3,t4;double kp,ep,exp(),delta;
CLET(v,*z);CTREAL(v,v,pi*.5/omega);
/*printf(" omega= %e m= %e eta=%e k=%e\n",omega,m,eta,k);*/
CMULT(dum1,*z,*z);CTREAL(dum1,dum1,eta/(2.*omega));
cexp(&dum1,&dum2);
/*printf(" before tek\n");*/
tek(0,1.-m,&kp,&ep);
/*printf(" after tek\n");*/
delta= g2*g2*g2-27.*g3*g3;
if(delta>0.)
	{CMPLX(q,exp(-pi*kp/k) ,0.);}
else
	{CMPLX(q,0.,exp(-pi*.5*kp/k));}
/*printf(" before ctheta v= %e %e q=%e %e\n",v.x,v.y,q.x,q.y);*/
ctheta(&v,&q,&t1,&t2,&t3,&t4);
/*printf(" after ctheta\n");*/
CMULT(v,t1,dum2);/*v= exp()*t1(v)*/
CMPLX(dum1,0.,0.);
/*printf(" before theta(0)\n");*/
ctheta(&dum1,&q,&t1,&t2,&t3,&t4);
/*printf(" after theta(0)\n");*/
CDIV(dum1,v,t2);CDIV(dum2,dum1,t3);CDIV((*ans) ,dum2,t4);
CTREAL (*ans,*ans,2.*omega/pi);
return 0;
}

zetaw(z,ans,k,omega,eta,m,g2,g3)
struct complex *z,*ans;double k,omega,eta,m,g2,g3;
{struct complex v,dum1,dum2,q,t1,t2,t3,t4;
double kp,ep,exp(),pow(),delta,mp,op,et,ehp,g2n,g3n;int n;
/* convergence trouble in cthld if |Im(z)| big */

/* try to use period relations to reduce Im z */
if(abs(z->y)>2.)
	{
	tek(0,1.-m,&kp,&ep);op=omega*kp/k;
	delta= g2*g2*g2-27.*g3*g3;
	if( delta>0.)	n=(int)((z->y)/(2.*op));
	else
		{/*op=omega2',want complex omega*/
		n=(int)(z->y/(op));
		}
/*printf(" reducing? op=%e n %d\n",op,n);*/
	if(n)
		{
		CLET(v,*z);
		if(delta>=0.)
			{v.y-=n*2.*op;CMPLX(dum2,0.,(eta*op-.5*pi)/omega);}
		else
			{v.y-=n*op;v.x-=n*omega;
			CMPLX(dum2,0.,(eta*op-pi)/omega);}
		zetaw(&v,&dum1,k,omega,eta,m,g2,g3);
		CTREAL(dum2,dum2, 2.*n);
		CADD(*ans,dum1,dum2);
		return 0;
		}
	}
if(abs(z->y)>.5 && omega>1.1)/* if omega not big no use only reduce once*/
	{/* for big |Im(z)| use */
	g2n=g2/pow(omega,4.);g3n=g3/pow(omega,6.);
/*printf(" reducing omega=%e,g2n g3n %e %e\n",omega,g2n,g3n);*/
	CLET(v,*z);CTREAL(v,v,1./omega);
	weier(&v,g2n,g3n,&t1,&t2,&mp,&kp,&ep,&op,&et,&ehp);
	zetaw(&v,&t3,kp,op,et,mp,g2n,g3n);
/*printf(" redu return\n");*/
	CTREAL(t3,t3,1./omega);
	CSET(ans,t3);
	return 0;
	}

CLET(v,*z);CTREAL(v,v,pi*.5/omega);
tek(0,1.-m,&kp,&ep);
delta= g2*g2*g2-27.*g3*g3;
printf(" g2, g3 delta %e %e %e\n",g2,g3,delta);
if(delta>0.)
	{CMPLX(q,exp(-pi*kp/k) ,0.);}
else
	{CMPLX(q,0.,exp(-pi*.5*kp/k));}
CSET(ans,*z);CTREAL(*ans,*ans,eta/omega);
ctld(&v,&q,&t1,&t2,&t3,&t4);
CTREAL(t1,t1,pi/(2.*omega));
CADD(*ans,*ans,t1);
return 0;
}

