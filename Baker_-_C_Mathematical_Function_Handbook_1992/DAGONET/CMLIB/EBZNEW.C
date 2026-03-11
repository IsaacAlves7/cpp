/*  Riemann Zeta (real arg.)  and related functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

zeta		Riemann Zeta
bernoulli(n)	Bernoulli numbers
bernpoly(n,x)	Bernoulli polynomial
euler		Euler numbers
debye		integral
dilog
clausen	clausen integral (Ch. 27)
*/

#include "cmlib.h"
#include "complex.h"
#include "protom.h"
#include <stdio.h>
#define itmax 1000
extern int itkt;

/* ibkeven is the integer power such that x^i is more cheaply
calculated using pow() than x*x*... i times*/
int ibkeven=11,dbkeven=3,lbkeven=4;/*default guess 387/386 TURBOC*/

double zeta(r) double r;
{
double pow(),bernoulli(),term,sum,g,sin(),gamma(),values[11]=
	{1.6449340668,1.20205690316,1.0823232337,1.0369277551,
	1.017343061984,1.0083492773819,1.00407735619,1.0020083928,
	1.000994575127,1.000494188604,1.000246086553};
int i,j,k,ir,intlim,rint,m;
if(r==1.)return errorcode;
if(r==0.)return -.5;
ir=(int)r;
rint= (r==(double)ir);
if(rint &&ir>1)
	{
	if(ir<13) return values[ir-2];
	/* if even,could either continue on to sum, or
	zeta(2m)= -1^(m+1) (2pi)^2m B[2m]/(2*(2m)!)*/
	}
if(rint &&ir<0)
	{
	if( !((-ir)%2 ))return 0.;/*even*/
	ir= (1-ir);
	return -bernoulli(ir)/(ir);
	}
if(r<0.)
	{
	return zeta(1.-r)*pow(2.*pi,r-1.)*2.*sin(pi*.5*r)*gamma(1.-r);
	}
if(r>1.)
	{
	intlim= pow( 32767.,1./r);/* max int we can safely raise to power*/
	sum=0;
	rint= rint && (ir<intlim);
	for(i=1;i<itmax;i++)
		{
		if(rint &&(i<intlim))
			{ m=k=(i<<1)-1 ;
			for(j=2;j<=ir;j++)m*=k;
			term=m;
			}
		else if(i<32767)term= pow( (double)((i<<1)-1),-r);
		/* if rint, ir<5 use long int possible*/
		else term= pow( 2.*((double)i)-1.,-r);
		sum+=term;
		}
	g=pow(2.,r);
	return sum*g/(g-1.);
	}
/* 0<r<1*/
g=-1.;sum=1.;
for(i=2;i< 20000;i++)
	{
	term= pow((double)i,-r);
	sum+=g*term;
	g=-g;
	if(abs(term)<1.e-5 || abs(term)<1.e-5*abs(sum))break;
	}
return sum/(1.-pow(2.,1.-r));
}

double bernoulli(n) int n;/* n even 0<=n*/
{
double r,sum,term,x,b0=1.,b1=-.5,b2=.166666666,g;
int i;
if(n==0)return b0;
if(n==2)return b2;
if(n==1)return b1;
if(n%2==1)return 0.;
r=(double) n;
if(n>3)
	{
	sum=0;
	for(i=1;i<itmax;i++)
		{
		term= pow( (double)(i) ,-r);
		sum+=term;
		if(abs(term)<5.e-7 || abs(term/sum)<5.e-7)break;
		}
	g=pow(2.*pi,r);
	x=1.;
	if((n>>1)%2 ==0)x=-1.;
	return 2.*sum*gamma(r+1.)/g*x;
	}
return 0.;
}

double bernpoly(n,x) double x;int n;
{int i; double sum,fact;
fact=1.;
sum= bernoulli(0);
for(i=1;i<=n;i++)
	{
	fact*=((double) (n+1-i))/((double)(i));
/*printf(" i %d sum %d fact %e\n",i,sum,fact);*/
	sum=sum*x+ fact*bernoulli(i);
	}
return sum;
}

double eulerpoly(n,x) int n; double x;
{int m;double y,z;
m=n+1;y=x*.5;
return 2./m*(bernpoly(m,x)-bernpoly(m,y)*pow(2.,(double)m));
}


double euler(n) int n;
{
double e[7]={1.,-1.,5.,-61.,1385.,-50521.,2702765.};
double r,g,sum,term,x;
int i;
itkt=0;
if(!n)return e[0];
if(n%2==1 ) return 0./*no such*/;
if( n<=12)
	{i=n>>1;return e[i];}
r=(double) (n);
if(n>3)
	{
	sum=0;
	x=1.;
	for(i=1;i<itmax;i++)
		{
		term= x*pow( (double)((i<<1)-1) ,-(r+1.));
		sum+=term;
		x=-x;
		itkt=i;
		if(abs(term)<1.e-10 || abs(term/sum)<1.e-8)break;
		}
	g=pow(2./pi,r+1.);
	x=-1.;
	if((n>>1)%2 ==0)x=1.;
	/*printf(" euler itkt=%d %e %e\n",itkt,term,sum);*/
	return x*sum*g*2.*gamma(r+1.);
	}
return 0.;
}

double debye(x,n) int n; double x;
{
int k,l;
double nn,np,nm,sum,term,interm,fk,factor,tt;
/*printf(" enter debye x=%le n=%d\n",x,n);*/
nn=n;np=nn+1.;
sum=0.;
factor=x*x;
term=factor;
if(x<.15 )
	{
	for(k=1;k<itmax;k++)
		{
		tt=term*bernoulli(2*k)/(2*k+n);
		sum+=tt;
		if(abs(tt)<1.e-5 || abs(tt/sum)<1.e-5)break;
		term*=factor/(2*(k+1)*(2*k+3));
		}
	return pow(x,nn)*(1./n-.5*x/np+sum);
	}

fk= exp(-x);
factor=fk;
for(k=1;k<itmax;k++)
	{
	term=0.;nm=nn;
	interm= pow(x, (double)n)/k;
	for(l=0;l<=n;l++)
		{
		term+=interm;
		interm*=(nm/(k*x));
		nm--;
		}
	term*= factor;
	factor*=fk;
	sum+=term;
	if(abs(term)<1.e-5 || abs(term/sum)<1.e-5)break;
	}
return zeta( np) *gamma(np)-sum;
}

double dilog(x) double x;
{
double log(),t;
if(x!=0.)t=log(x);
else t= -1.e60;
/* x=0 -> t= -infinity. do not use errorcode here, as user
might change it from -1.e60 */

if(x>1.)
	{
	return -t*t*.5-debye(t,1);
	}
return debye(-t,1);
}


double clausen(theta) double theta;
{
int k;
double sum,term,sq,ans1,ans2,tt,t;
/*printf(" entered clausen theta= %le %le\n",theta, theta*180./pi);*/
if(theta==0.)return 0.;
if(theta<0. || theta>pi)return errorcode;
if(theta<=.5* pi )
	{
	sum=0.;
	sq=theta*theta;
	term=.25*theta*sq;
	for(k=1;k<itmax;k++)
		{
		tt=term*bernoulli(2*k)/((double)(k*((k<<1)+1)));
		sum+=tt;
/*printf(" debug sum term %e %e\n",sum,term);*/
		term*= -sq/((double)(2*(k+1)*(2*k+1)) );
		if((abs(tt)<5.e-7 || abs(tt)<1.e-6*abs(sum))  )break;
		itkt=k;
		}
	ans1=theta-theta*log(theta)+sum; return ans1;
	}
/* theta >= pi/2*/
/*printf(" debug > pi/2\n");*/

t= pi -theta;
sum=0.;
sq=t*t;
term=.25*t*sq;
for(k=1;k<itmax;k++)
	{
	tt=term*bernoulli(2*k)*(pow(2.,2.*k)-1.)/((double)(k*(2*k+1)));
	sum+=tt;
/*printf(" debug sum term %e %e\n",sum,tt);*/
	term*= -sq/(2.*(k+1)*(2*k+1) );
	if(abs(tt)<5.e-7 || abs(tt)<1.e-6*abs(sum))break;
	itkt=k;
	}
ans2=t*log(2.)-sum;
/*
sum=0.;
for(k=1;k<itmax;k++)
	{
	term= sin(k*theta)/((double)(k*k));
	sum+=term;
	printf(" debug sin sum term %e %e\n",sum,term);
	if((abs(term)<5.e-7 || abs(term)<1.e-6*abs(sum))&& k>40)break;
	itkt=k;
	}
ans1=sum;*/
/*printf(" clausen sin sum=%le %le\n",ans1,ans2);*/
return ans2;
/*return ans1;*/
/* leads to infinite loop theta>pi/2
t=theta*.5;
return 2.*(clausen(t*.5)-clausen(pi-t));
*/
}


double zeta1(z) double z;
{double bigphi();/* for faster see ebz.c*/
return bigphi(1.,z,1.);
}


double zeta2(z,a)double z,a;
{/* "Hurwitz" zeta function*/
int i,top,topl=12,inta;/*topl chosen for efficiency*/
double zeta(),sum,pow(),term,tol=1.e-6,denom;
if(a==0.)return errorcode;
if(z==0.)return .5-a;
if (((z-(double)((int)z))==0.) && z<0.)
	{/* z negative integer*/
	return -bernpoly((int)(1.-z),a)/(z+1.);
	}
inta=((a-(double)((int)a))==0.);
if(  inta && a>0.)
	{/* a is a positive  integer*/
	if(a==1.)return zeta(z);
	top= (int)a;
	if(top<topl)
		{
		sum=0.;
		for(i=1;i<top;i++)
			{
			term=pow(((double)i) ,-z);
			sum+=term;
			if( abs(term)<tol*abs(sum))break;
			}
/*		printf(" answer shortcut %e\n",zeta(z)-sum);*/
		return zeta(z)-sum;
		}
	}
sum=0.;
for(i=0;i<itmax;i++)
	{
	denom=a+((double)i);
	if(denom!=0.)term=pow( denom,-z);
	sum+=term;
	if( abs(term)<tol*abs(sum))break;
	}
return sum;
}


double bigphi(z,s,v) double z,s,v;
{
double sum,term,pow(),denom,power,tol=1.e-6;
int i;
sum=0.;
power=1.;
for(i=0;i<itmax;i++)
	{
	denom=v+((double)i);
	term=pow(denom ,-s)*power;
	sum+=term;
	if( abs(term)<tol*abs(sum))break;
	power*=z;
	}
return sum;
}

double betacat(s) double s;
{double pow();/* Catalan beta= sum k=0 to inf. (-1)^k (2k+1)^-s*/
if(s==1.) return pi*.25;
if(s==2.) return .915965594177219;
if(s==3.)return pi*pi*pi*.03125;
return bigphi(-1.,s,.5)*pow(2.,-s);}

double lambda(s) double s;
{double pow(),zeta();/* sum k=0 to infinity (2k+1)^-s */
return (1.-pow(2.,-s))*zeta(s);
}

double eta(s) double s;
{double pow(),zeta();/* sum from k=1 to infinity -1^k k^-n */
return zeta(s)*(1.-pow(2.,1.-s));
}

double ifermi(mu,s)double mu,s;
{/* integral from 0 to infinity of k^s/(exp(k-mu)+1) dk */
double exp(),gamma(),e;e=exp(mu);
return e*gamma(1.+s)*bigphi(-e,1.+s,1.);
}

