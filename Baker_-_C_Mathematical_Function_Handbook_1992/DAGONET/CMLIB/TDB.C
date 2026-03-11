/*  
test driver for bessel functions
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

extern double c1a,c2a;
extern double sphb[200];
extern int itkt;

FILE *intk,*inti,*intj,*intkt,*intit,*inty,*intjt,*ints,*inta,*aff;

main()
{int i,m,n,rrr;double x,y,t,q,p,nu,xx,rr,r;
int check;
struct complex z,j,yy,h2,jp,yp,h2p;

intk=fopen("intk.out","w");fprintf(intk," 2\n");
ints=fopen("ints.out","w");  fprintf(ints," 3\n");
intit=fopen("intit.out","w");  fprintf(intit," 2\n");
inti=fopen("inti.out","w");  fprintf(inti," 2\n");
intkt=fopen("intkt.out","w");fprintf(intkt," 2\n");
intjt=fopen("intjt.out","w");fprintf(intjt," 2\n");
intj=fopen("intj.out","w");fprintf(intj," 2\n");
inty=fopen("inty.out","w");fprintf(inty," 2\n");
inta=fopen("inta.out","w");fprintf(inta," 2\n");
aff=fopen("af.out","w");fprintf(aff," 5\n");
for(n=0;n<10;n++)
	{
	x=(n+1)*.1;
	p=ybes(x,1.);q=y1(x);
printf(" Ybes x %e y=%e %e\n",x,p,q);
	}
for(n=0;n<10;n++)
	{
	x=(n+1)*.1;
	p=ibes(x,1.);
printf(" Ibes x %e i=%e \n",x,p);
	}
for(n=0;n<10;n++)
	{
	x=(n+1);
	p=ibes(x,1.);
printf(" Ibes x %e i=%e \n",x,p);
	}
for(n=0;n<10;n++)
	{
	x=(n+1)*.1;
	p=kbes(x,.5);
	t=ibes(x,.5);
printf(" i=%e",t);
	q=ybes(x,.5);
printf(" x %e k,y,i=%e %e %e\n",x,p,q,t);
	}
for(n=0;n<10;n++)
	{
	x=(n+1);
	p=kbes(x,.33333);
printf("  k(%e,1/3)=%e \n",x,p);
	}
for(n=0;n<10;n++)
	{
	x=(n+5);
	p=jbes(x,1.);
	q=j1(x);
printf(" x %e j1=%e %e \n",x,p,q);
	}
/* tables of bessel functions*/
for(n=0;n<10;n++)
	{
	x=.1*(n+1);
	p=jt(x);
	nu=jbes(x,.333333);
printf(" x %e jt=%e j=%e\n",x,p,nu);
	}

for(n=0;n<10;n++)
	{
	x=(n+1);
	p=jbes(x,0.);
	nu=jbes(x,1.);
	t=jbes(x,3.);
	q=jbes(x,4.);
	printf(" x %e j0=%e j1=%e \n  j3=%e j4=%e\n",x,p,nu,t,q);
	}
printf(" jas,j(2.,1/3)=%e %e\n",jas(2.,.333333),jbes(2.,.333333333));
printf(" jas,j(3.,1/3)=%e %e\n",jas(3.,.333333),jbes(3.,.333333333));
printf(" jas,j(4.,1/3)=%e %e\n",jas(4.,.333333),jbes(4.,.333333333));
printf(" jas,j(5.,1/3)=%e %e\n",jas(5.,.333333),jbes(5.,.333333333));
printf(" jas,j(10.,1/3)=%e %e\n",jas(10.,.333333),jbes(10.,.333333333));
printf(" jas,j(15.,1/3)=%e %e\n",jas(15.,.333333),jbes(15.,.333333333));
printf(" jas,j(20.,1/3)=%e %e\n",jas(20.,.333333),jbes(20.,.333333333));
printf(" jas,j(2.,0.)=%e %e\n",jas(2.,0.),jbes(2.,0.));
printf(" jas,j(3.,0)=%e %e\n",jas(3.,0.),jbes(3.,.0));
printf(" jas,j(4.,0)=%e %e\n",jas(4.,0.),jbes(4.,.0));
printf(" jas,j(5.,0)=%e %e\n",jas(5.,0.),jbes(5.,.0));
printf(" jas,j(10.,0)=%e %e\n",jas(10.,0.),jbes(10.,.0));
printf(" jas,j(15.,0)=%e %e\n",jas(15.,0.),jbes(15.,.0));
printf(" jas,j(20.,0)=%e %e\n",jas(20.,0.),jbes(20.,.0));
/* functions based upn complex bessel functions*/
while(1)
	{
	printf(" enter n z-complex\n");scanf("%d%le%le",&n,&q,&r);
	if(q==0. && r==0. && n==0)break;
	z.x=q;z.y=r;
	bessel(n,&z,&j,&yy,&h2,&jp,&yp,&h2p,&check);
	printf(" J,Y=");printc(&j);printc(&yy);printf("\n");
	printf(" J',Y'=");printc(&jp);printc(&yp);printf("\n");
	printf(" H2,H2'=");printc(&h2);printc(&h2p);printf("\n");
	ibess(&z,&j,n); printf(" Ibess=");printc(&j);printf("\n");
	kbess(&z,&j,n); printf(" Kbess=");printc(&j);printf("\n");
	if( r==0.)
		{printf(" jn=%le %le\n", jint(n,q),jn(q,n));
		kelvin(n,(double)q,&j,&yy);
		printf(" kelvin be=");printc(&j);
		printf(" kelvin ke=");printc(&yy);
		printf("\n");
		}
	}

printf("\n\n Chapter 10: Spherical Bessel functions \n\n");


/* spherical bessel functions:*/
printf(" Spherical bessel functions:\n");
while(1){
printf(" spherical i:enter n, x\n");
scanf("%d%le",&m,&x);
if(x<0.)break;
printf(" %e\n",sinb(x,m));
for(i=0;i<=abs(m);i++)
	printf(" i[%d]=%e\n",i,sphb[i]);
if(x>0.)
	{
	if(m>0)printf(" ck: i[0]=%le\n", sinh(x)/x);
	if(m>1)printf(" ck: i[1]=%le\n", (cosh(x)-sinh(x)/x)/x);
	if(m>2)printf(" ck: i[2]=%le\n",
	  (-3.*cosh(x)/x+sinh(x)*(1.+3./(x*x)))/x);
	if(m<=-1)printf(" ck: i(-1) pro I[-1/2]=%le\n",
		cosh(x)/x);
	if(m<=-2)printf(" ck: i(-2) pro I[-3/2]=%le\n",
		(cosh(x)*(-1./x)+sinh(x))/x);
	if(m<=-3)printf(" ck: i(-3) pro I[-5/2]=%le\n",
		(cosh(x)*(1.+3./(x*x))-3./(x)*sinh(x))/x);
	}
}
while(1){
printf(" k:enter n, x\n");
scanf("%d%le",&m,&x);
if(x<0.)break;
printf(" %e\n",skn(x,m));
for(i=0;i<=m;i++)
	printf(" %d %e\n",i,sphb[i]);
}
while(1){
printf(" j:enter n, x\n");
scanf("%d%le",&m,&x);
if(x<0.)break;
printf(" %e\n",sjn(x,m));
for(i=0;i<=m;i++)
	printf(" %d %e\n",i,sphb[i]);
}
while(1){
printf(" y:enter n, x\n");
scanf("%d%le",&m,&x);
if(x<0.)break;
printf(" %e\n",syn(x,m));
for(i=0;i<=m;i++)
	printf(" %d %e\n",i,sphb[i]);
}

printf("\n\n Chapter 10: Airy Functions & Integrals\n\n");

for(n=0;n<10;n++)
	{
	x=.1*(n+1);
	p=ai(x);
	nu=bi(x);t=ai(-x);q=bi(-x);
printf(" x %e Ai,Bi=%e %e %e %e\n",x,p,nu,t,q);
	}
for(n=0;n<10;n++)
	{
	x=(n+1);
	p=ai(x);
	nu=bi(x);t=ai(-x);q=bi(-x);
printf(" x %e Ai,Bi=%e %e %e %e\n",x,p,nu,t,q);
	}
for(n=0;n<10;n++)
	{
	x=(n+1);
	p=dai(x);nu=dbi(x);
	t=dai(-x);q=dbi(-x);
printf(" x %e Ai',Bi'=%e %e %e %e\n",x,p,nu,t,q);
	}

for(n=0;n<10;n++)
	{
	x=.1*(n+1);
	p=ai(x);
	nu=bi(x);t= c1a*smallf(x)-c2a*smallg(x);
	q=sqrt(3.)*(c1a*smallf(x)+c2a*smallg(x));
printf(" x=%e Ai=%e Ai(via f,g)=%e \n   Bi=%e via(f,g)=%e\n"
	,x,p,t,nu,q);
	}
for(n=0;n<100;n++)
	{
	x=.1*(n+1)-5.;
	p=ai(x);
	nu=bi(x);
	fprintf(aff," %le %le %le\n",x,p,nu);
	}
for(n=0;n<10;n++)
	{
	x=(n+1);
	p=IAi(x);nu=IAi(-x);t=IBi(x);q=IBi(-x);
	printf(" x %f Airy integrals=%e %e %e %e\n",x,p,nu,t,q);
	}
for(n=0;n<50;n++)
	{
	x=(n+1)*.2;
	p=IAi(x);nu=IAi(-x);t=IBi(x);q=IBi(-x);
	fprintf(inta," %f %e %e %e %e\n",x,p,nu,t,q);
	}

printf("\n\n Chapter 11: Integrals\n\n");

	infinite_loop
		{
		printf(" enter x>=0, r for bickley x==0 analytic \n");
		scanf("%le%le",&xx,&rr);
		if(xx<0.)break;
		printf(" x, r %le %le\n",xx,rr);
		printf("Ki=%le\n",bickley(xx,rr));
		if(xx==0.)
			{rrr=rr;
			m=rrr%2;
			if(m) {rr=(rr-1.)*.5;
					q= pi*.5*gamma(rr+.5)/(gamma(.5)*gamma(rr+1.));}
			else {rr*=.5;q=gamma(rr)*gamma(1.5)/gamma(rr+.5);}
			printf(" analytic %le\n",q);
			}
		}

infinite_loop
		{
		printf(" enter x, r,n[int] for Jrn");scanf("%le%le%d",&xx,&rr,&m);
		if(xx==0. && rr==0. && m==0)break;
		printf(" x, r %le %le n=%d\n",xx,rr,m);
		printf(" rth integral:Jn=%le\n",Jrn(xx,rr,m));
		}

/* integrals of bessel functions:*/
printf(" integrals of bessel functions\n");
for(i=1;i<20;i++)
	{x=i*.5;
	y=0.;if(x>=5.)y=ii0tas(x);
	t=ii0m1t(x);printf(" I0/t %e %e %e %e %e %d\n",x,t,t*exp(-x),y,y*exp(-x),itkt);
	}
for(i=1;i<40;i++)
	{x=i*.25;
	t=ii0m1t(x);
	fprintf(intit," %le %le\n",x,t);
	}
for(i=1;i<10;i++)
	{x=i*.5;t=iy0t(x);printf("int Y0/t %e %e \n",x,t);}
for(i=1;i<50;i++)
	{x=i*.1;t=iy0t(x);
	fprintf(inty,"%e %e \n",x,t);
	}
for(i=0;i<10;i++)
	{x=i*.5;t=ki(x);printf(" int K0 %e %e \n",x,t);
	}
for(i=0;i<50;i++)
	{x=i*.1;t=ki(x);printf(" int K0 %e %e \n",x,t);
	fprintf(intk,"%e %e \n",x,t);
	}

for(i=1;i<20;i++)
	{x=i*.5;t=ik0t(x);
	y=0.;if(x>=4.)y=ik0tas(x);
	printf("int K0/t %e %e %e %e %e %d\n",x,t,t*x*exp(x),y,y*x*exp(x),itkt);
	fprintf(intkt," %le %le\n",x,t);
	}

	for(i=1;i<10;i++)
	{x=i*.5;t=ij0t(x);printf("int J0/t %e %e \n",x,t);
	}
	for(i=1;i<50;i++)
	{x=i*.5;t=ij0t(x);
	fprintf(intjt," %le %le\n",x,t);
	}

for(i=1;i<25;i++)
	{x=i*.5;t=ji(x);nu=jin(x,1);printf(" int J0 %le %le %le\n",x,t,nu);
	fprintf(intj," %le %le\n",x,t);
	}




infinite_loop
	{printf(" enter x for integral i0i y0i\n");
	scanf("%le",&t);if(t==0.)break;
	printf(" integrals I0 %le Y0 %le\n",i0i(t),y0i(t));
	}

for(i=0;i<50;i++)
	{t=(i+1)*.2;
	fprintf(inti,"%le %le %le\n",t,i0i(t),y0i(t));
	}


printf("\n\n Chapter 12: Struve and Anger-Weber functions & Integrals\n\n");

infinite_loop
		{printf(" enter x, nu\n");scanf("%le%le",&x,&nu);
		if(x==0. && nu==0.)break;
		printf(" H= %le L=%le\n", StruveH(nu,x),StruveL(nu,x));
		}

for(i=0; i<50;i++)
	{t=i*.1;
	rr=i0i(t)-ModStruveI(t);
	if(t>5.)rr-= 2./pi*log(t);
	printf("x=%le int H %le H/t*2/pi= %le, \n  int L %le f2=%le\n"
		,t,StruveI(t),StruveIot(t)*2./pi,ModStruveI(t),rr);
	fprintf(ints,"%le %le %le\n",t,StruveI(t),StruveIot(t));
	}


printf(" anger weber\n");
aw(.5,.5,&p,&t);
printf(" nu=x=.5, %e %e\n",p,t);
printf(" anger weber\n");
aw(.5,10.,&p,&t);
printf(" nu=.5 x=10, %e %e\n",p,t);

t=0.;
for(n=0;n<10;n++)
	{
	x=.2*(n+1);aw(t,x,&p,&nu);
	q=-StruveH(0.,x);
	printf(" x %e nu=%e J=%le j0=%le \n   E=%le %le\n"
		,x,t,p,j0(x),nu,q);
	}
t=1.;
for(n=0;n<10;n++)
	{
	x=.2*(n+1);aw(t,x,&p,&nu);
	q= 2./(pi)-StruveH(1.,x);
	printf(" x %e nu=%e J=%le j1=%le \n   E=%le %le\n"
		,x,t,p,j1(x),nu,q);
	}
for(n=0;n<10;n++)
	{ /* ck asymptotics*/
	x=(n+1);aw(t,x,&p,&nu);
	q= 2./(pi)-StruveH(1.,x);
	printf(" x %e nu=%e J=%le j1=%le \n E=%le %le\n"
		,x,t,p,j1(x),nu,q);
	}
infinite_loop
	{
	printf(" enter int m, complex z for integral Anger Weber\n");
	scanf("%d%le%le",&n,&q,&p);
	if(q==0. && p==0. && n==0.)break;
	yy.x=q;yy.y=p;
	iaw(n, &yy,&z);
	printf(" answer= %le %le\n", z.x,z.y);
	}

return 0;
}
