/*  test driver for Legendre functions of chapter 8.

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
#include <stdio.h>
#include <alloc.h>
#include "cmlib.h"
#include "protom.h"

extern double leg[100];

main()
{
double s2,f21(),q0,q1;
float   x,tau; double z,y,b,p[200]; int i,l,n,m,d=7,real,rtcode;
double q[200],r[200];

while(1)
	{
	printf(" P: enter nmax, m, x\n");
	scanf("%d%d%le",&n,&m,&z);
	if(n==0 && m==0 && z==0.)break;
	legptable(z,n,m);
	for(i=0;i<=n;i++)printf(" %d %le\n",i,leg[i]);
	}

/* exercise Gaut.c and added routines*/
while(1)
	{
	printf(" P: enter n, mmax, x>1 \n");
	scanf("%d%d%e",&n,&m,&x);z=x;
	if(n==0 && m==0 && x==0.)break;
	rtcode=leg1(z,n,m,p);
	if(rtcode)continue;
	for(i=0;i<=m;i++)printf(" %d %le\n",i,p[i]);
	}
while(1)
	{
	printf(" Q:enter nmax, m, x\n");
	scanf("%d%d%e",&n,&m,&x);z=x;
	if(n==0 && m==0 && x==0.)break;
	rtcode=leg2(z,m,n,d,p);
	if(rtcode)continue;
	for(i=0;i<=n;i++)printf(" %d %le\n",i,p[i]);
	}
while(1)
	{
	printf(" Q:enter n, mmax, x\n");
	scanf("%d%d%e",&n,&m,&x);z=x;
	if(n==0 &&m==0 && x==0.)break;
	rtcode=leg3(z,n,m,d,p);
	if(rtcode)continue;
	for(i=0;i<=m;i++)printf(" %d %le\n",i,p[i]);
	}
while(1)
	{
	printf(" P:enter alpha mmax x\n");
	scanf("%e%d%e",&tau,&m,&x);z=x;y=tau;
	if(tau==0. &&m==0 && x==0.)break;
	rtcode=legend1(z,y,m,d,p);
	if(rtcode)continue;
	for(i=0;i<=m;i++)printf(" %d %le\n",i,p[i]);
	}

while(1)
	{
	printf(" P:enter alpha+nmax m x\n");
	scanf("%e%d%d%e",&tau,&n,&m,&x);z=x;y=tau;
	if(tau==0. && n==0 &&m==0 && x==0.)break;
	rtcode=legend2(z,y,m,n,d,p);
	if(rtcode)continue;
	for(i=0;i<=n;i++)printf(" %d %le\n",i,p[i]);
	}

printf(" toroidal chk Q:\n");
for(i=0;i<10;i++)
	{
	z=1.+(i+1)*.1;
	toroidal(z,0,1,d,p);
	/*y= sqrt((z-1.)/(z+1.));*/
	s2=(2./(z+1.));            y=sqrt(s2);
	tek(0,s2,&q0,&q1);
	q1= z*y*q0-sqrt(2.*(1.+z))*q1;
	printf(" z %f Q-1/2,Q1/2=%le %le; elliptic I: %le %le \n"
		,z,p[0],p[1],y*q0,q1);
	}

printf(" toroidal chk\n");
for(i=0;i<10;i++)
	{
	b=1.+(i+1)*.1;
	 z= cosh(b);
	toroidal(z,0,1,d,p);
	/*  7.10.1 of Lebedev*/
	q0= sqrt(pi)*gamma(.5)/gamma(1.)*exp(-.5*b)*f21(.5,.5,1.,exp(-2.*b));
	q1= sqrt(pi)*gamma(1.5)/gamma(2.)*exp(-1.5*b)*f21(1.5,.5,2.,exp(-2.*b));

	printf(" z %f %le %le Q-1/2,P1/2 ck: %le %le\n",z,p[0],p[1],q0,q1);
	}

while(1)
	{
	printf(" Q toroidal:enter m nmax x\n");
	scanf("%d%d%e",&m,&n,&x);z=x;y=tau;
	if(m==0 &&n==0 && x==0.)break;
	toroidal(z,m,n,d,p);
	for(i=0;i<=n;i++)printf(" %d %le\n",i,p[i]);
	}
while(1)
	{
	printf(" P toroidal:enter n x\n");
	scanf("%d%e",&n,&x);z=x;
	if(n==0 && x==0.)break;
	printf(" toroidal P=%le\n",Ptoroidal(n,z));
	}

while(1)
	{printf(" enter 1 for  Q(x) 0 Q(ix) -1 stop\n");scanf("%d",&real);
	if(real<0 ) break;
	printf(" enter nu,x for Q\n");scanf("%le%le",&b,&z);
	printf(" x=%f nu=%f %le\n",z,b,qnu(z,b,real));
	}

while(1)
	{
	printf("  conical enter tau, theta\n");
	scanf("%e%e",&tau,&x); z=x;b=tau;
	if(z==0. && b==0.)break;
	printf(" ans=%le theta=%le tau=%le\n",conicalt(z,b),z,b);
	}

printf(" asympt. tau conical ck\n");
for(i=1;i<10;i++)
	{
	/* small ztheta= pi*i/10.;
	z=cos(theta);
	rtcode=conical(z,tau,0,d,p);
	if(rtcode)continue;
	s2= sin(theta*.5);s2*=s2;
	y= 1. +(4.*tau*tau+1.)/(4.)*s2
	 + (4.*tau*tau+1.)*(4.*tau*tau+9.)/(4.*16.)*s2*s2;
	printf(" ans=%le vs %le\n",p[0],y);*/
	b= 2.;
	z=cosh(b);
	tau =10.*i;
	rtcode=conical(z,tau,0,d,p);
	if(rtcode)continue;
	q0=sqrt(2./(pi*tau*sinh(b)))*sin(tau*b+.25*pi);
	printf(" z tau %f %f ans %le %le\n",z,tau,p[0],q0);
	}

while(1)
	{
	printf(" conical:enter tau nmax x\n");
	scanf("%e%d%e",&tau,&n,&x);z=x;y=tau;
	if(tau==0. &&n==0 && x==0.)break;
	rtcode=conical(z,y,n,d,p);
	if(rtcode)continue;
	for(i=0;i<=n;i++)printf(" %d %le\n",i,p[i]);
	}

/* leg.c: Hernden*/
while(1)
	{
	printf(" enter l,m,r=1 real(<0 to quit),x\n");
	scanf("%d%d%d%f",&l,&m,&real,&x);
	if(real<0)break; z=x;
	printf(" %le\n",legendrea(m,l,z,real));
	}
while(1)
	{
	printf(" enter l,m,r=1 real,x>0. real<0 to quit\n");
	scanf("%d%d%d%f",&l,&m,&real,&x);
	if(real<0)break; z=x;   qleg(m,l,z,real,r,q);
	for(i=0;i<=l;i++)printf(" %d r=%le q=%le\n",i,r[i],q[i]);
	}
/* leg.c: other */
while(1)
	{
	printf(" plm enter l,m,x real,x<1000.\n");
	scanf("%d%d%le",&l,&m,&z);
	if(z>1000.)break; 
	printf(" for real x: %le for imag x, m=0: %le\n",plm(l,m,z),pli(l,z));
	}
while(1)
	{
	printf(" qlm enter l,m,x real,x<1000.\n");
	scanf("%d%d%le",&l,&m,&z);
	if(z>1000.)break; 
	printf(" for real x: %le for imag x, m=0: %le\n",qlm(l,m,z),qli(l,z));
	}
while(1)
	{
	printf(" pmu,nu & qmu,nu: enter mu,nu,x real,x<1000.\n");
	scanf("%le%le%le",&y,&b,&z);
	if(z>1000.)break; 
	printf(" P,Q: %le %le\n",pmunu(y,b,z),qmunu(y,b,z));
	}

/* Mehler*/
while(1)
	{
	printf(" enter x,z,mm Mehler x>1000 to quit\n");
	scanf("%le%le%d",&z,&y,&m);
	if(z>1000.)break;
	printf(" Mehler=%le for z=0: %le\n",Mehler(z,y,m),Mehler0(z,y));
	}

}