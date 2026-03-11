/*  Elliptic Integral of Third Kind

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

double e3(n,phi,m,sig) double n,phi,m,sig;
{
int choice;
double eps=1.e-6,sin(),cos(),pye,ms,cs,new,p1,y,dp,x,dnr,ddr,t1,t2,t3,
	atan(),tan(),log(),sqrt(),fp,ep,sinp,secp,tanp,snp,snm,sqn,sqm,tpd,tpn,
	ca,tp,rhs,py,dpn,dpd,s2p,sgn;
choice=0;
if(n==0.)
	{if(m==0.)return phi;tef(phi,m,sig,&fp,&ep);return fp;
	}
if(n>0. && n<m)choice=1;
if(n>1.)choice=2;
if(n>m && n<1.)choice=3;
if(n<0.)choice=4;
if( abs(n)<=eps && m>eps)choice=5;
if( abs(n)<=eps && m<=eps)choice=6;
if(m<eps && abs(n)>eps)choice=7;
if(m==1. && n!=1.)choice=8;
ms=sqrt(m);
if( abs(n-ms)<=eps || abs(n+ms)<=eps)choice=9;
cs=sqrt(1.-m);
if( abs(n-1.-cs)<eps || abs(n-1.+cs)<eps)choice=11;
if(abs(n-m)<=eps)choice=13;

switch (choice)
{
case 1:
	e3nlm(n,phi,m,sig,&pye);
	return pye;
case 2:
	new=m/n;
	p1=sqrt((m-1.)*(1.-new));
	tef(phi,m,sig,&fp,&ep);
	y=sin(phi);
	dp=sqrt(1.-m*y*y);
	x=p1*tan(phi);
	dnr=dp+x;
	ddr=dp-x;
	t3=.5/p1*log(dnr/ddr);
	e3nlm(new,phi,m,sig,&pye);
	pye= -pye+fp+t3;
	return pye;
case 3:
	e3mln(n,phi,m,sig,&pye);
	return pye;
case 4:
	new=(m-n)/(1.-n);
	p1=sqrt(-n*new);
	e3mln(new,phi,m,sig,&pye);
	t1=sqrt((1.-new)*(1.-m/new))*pye;
	tef(phi,m,sig,&fp,&ep);
	t2=m/p1*fp;
	x=sin(phi);
	dp=sqrt(1.-m*x*x);
	t3=atan(.5*p1*sin(2.*phi)/dp);
	pye=(t1+t2+t3)/sqrt((1.-n)*(1.-m/n));
	return pye;
case 5:
	tef(phi,m,sig,&pye,&ep);
	return pye;
case 6:	return phi;
case 7:
	if(n==1.)
		return tan(phi);
	if(n<1.)
		p1=sqrt(1.-n);
	else
		p1=sqrt(n-1.);
	return atan(p1*tan(phi))/p1;	
case 8:
	secp=1./cos(phi);
	sinp=sin(phi);
	tanp=sinp*secp;
	t1= log(tanp+secp);
	sqn=sqrt(n);
	x=sqn*sinp;
	snp=1.+x;
	snm=1.-x;
	t2=.5*sqn*log(snp/snm);
	return (t1-t2)/(1.-n);
case 9:
	tp=tan(phi);
	x=sin(phi);
	dp=sqrt(1.-m*x*x);
	sqm=sqrt(m);
	if(abs(n-ms)<=eps) sgn=-1.;
	else sgn=1.;
	y=1.+sgn*sqm;
	rhs=atan(y*tp/dp);
	py=rhs/y;
	tef(phi,m,sig,&fp,&ep);
	return .5*(py+fp);
case 11:
	tp=tan(phi);
	x=sin(phi);
	ca=sqrt(1.-m);
	dp=sqrt(1.-m*x*x);
	y=tp*dp;
	tpn=1.+y;
	tpd=1.-y;
	if(abs(n-1.+cs)<=eps)sgn=-1.;
	else sgn=1.;
	y=ca*tp;
	dpn=dp+y;
	dpd=dp-y;
	tef(phi,m,sig,&fp,&ep);
	t1= sgn*.5*log(tpn/tpd);
	t2=.5*log(dpn/dpd);
	t3=-sgn*(1.-sgn*ca)*fp;
	return (t1+t2+t3)/(2.*ca);
case 13:
	tef(phi,m,sig,&fp,&ep);
	p1=1./(1.-m);
	t2=m*p1;
	s2p=sin(2.*phi);
	y=sin(phi);
	dp=sqrt(1.-m*y*y);
	t1=p1*ep;
	t2 *=s2p/(2.*dp);
	return t1-t2;
/*if( abs(n-1.)<eps && m>eps)
	{
	}*/
default:
tef(phi,m,sig,&fp,&ep);
tp=tan(phi);
y=sin(phi);
p1=1./(1.-m);
dp=sqrt(1.-m*y*y);
pye=fp-p1*ep+p1*tp*dp;
return pye;
}
}

e3nlm(n,phi,m,sig,pye) double n,phi,m,sig,*pye;
{
int i;
double x,k,e,b,fe,ee,asin(),sqrt(),exp(),tan(),sin(),tol=1.e-7,
	xq,q,q2,m1,y,z,d1,v,ris,term,sum1,sum2,cotb,k1,e1,fp,ep,arg;
x=n/m;
e=asin(sqrt(x));
tef(e,m,sig,&fe,&ee);
tek(0,m,&k,&e);
b=pi*.5*fe/k;
cotb=1./tan(b);
m1=1.-m;
tek(0,m1,&k1,&e1);
tef(phi,m,sig,&fp,&ep);
arg=-pi*k1/k;
q=exp(arg);
v=pi*.5*fp/k;
xq=1.;
d1=sqrt(n/((1.-n)*(m-n)));
sum1=sum2=0.;
for(i=1;i<200;i++)
	{
	ris=i;
	xq*=q;
	if(xq<1.e-15)break;
	term=2.*xq*sin(2.*ris*v)*sin(2.*ris*b)/(ris*(1.-xq*xq));
	sum1+=term;
	if(i>1 && (abs(term)<tol)||( abs(term)/sum1 <sig)) break;
	}
y=sin(2.*b);
z=cos(2.*b);
xq=1.;
q2=q*q;
for(i=0;i<200;i++)
	{
	xq*=q2;
	if(xq<1.e-14)break;
	term=4.*xq*y/(1.+xq*(xq-2.*z)); 
	sum2+=term;
	if(i && (abs(term)<tol)||( abs(term)/sum2 <sig)) break;
	}
*pye=d1*(v*(cotb+sum2)-sum1);
return 0;
}

e3mln(n,phi,m,sig,pye) double n,phi,m,sig,*pye;
{
int i,iq;
double cosh(),sinh(),pow(),exp(),asin(),sqrt(),tanh(),atan(),
e,m1,fe1,e1,fp,ep,km1,em1,b,q,xx,term,
	ris,s2sv,sb,sb2,km,em,qs,arg,
	sum1,sum2,sum3,sgn,q2,q2s,d2,v,tol=1.e-7,qs2,mu,lm,thb,tv;
e=asin(sqrt((1.-n)/(1.-m)));
m1=1.-m;
tef(e,m1,sig,&fe1,&e1);
tek(0,m,&km,&em);
tef(phi,m,sig,&fp,&ep);
tek(0,m1,&km1,&em1);
b=.5*pi*fe1/km;
arg=-pi*km1/km;
q=exp(arg);
xx=fp/km;
if(xx>.99999)xx=.99999;
v=.5*pi*xx;
d2=sqrt(n/((1.-n)*(n-m)));
q2=q*q;
q2s=1.;
sum1=sum2=sum3=0.;
sgn=-1.;
for(i=1;i<200;i++)
	{
	ris=i;
	sgn=-sgn;
	q2s*=q2;
	s2sv=sin(2.*ris*v);
	sb2=2.*ris*b;
	if(abs(sb2)>10.)
		{
		sb=1.e10;
		if(sb2<-10.) sb=0.;
		}
	else sb=sinh(sb2);
	term=2.*sgn*q2s*s2sv*sb/(ris*(1.-q2s));
	sum1+=term;
	if(i>1 && (abs(term)<tol)||( abs(term)/sum1 <sig)) break;
	}
qs2=q;
for(i=1;i<200;i++)
	{
	ris=i;
	sb2=2.*ris*b;
	term=ris*qs2*sinh(sb2);
	sum2+=term;
	iq=(i<<1)+1;
	qs2*= pow(q,(double)iq);
	if(qs2<1.e-8)break;
	if(i>1 && (abs(term)<tol)||( abs(term)/sum2 <sig)) break;
	}
qs=q;
for(i=1;i<200;i++)
	{
	ris=i;
	sb2=2.*ris*b;
	term=2.*cosh(sb2)*qs;
	sum3+=term;
	iq=(i<<1)+1;
	qs2*= pow(q,(double)iq);
	if(qs2<1.e-8)break;
	if(i>1 && (abs(term)<tol)||( abs(term)/sum3 <sig)) break;
	}
thb=tanh(b);/* b as in my FTN or of sb2 math notes?*/
tv=tan(v);
lm=atan(thb*tv)+sum1;
mu=sum2/(1.+sum3);
*pye=d2*(lm-4.*mu*v);
return 0;
}