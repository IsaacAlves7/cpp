/*  test driver jacobian elliptic function and relatives

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"


extern double jtheta, amplitude,kimag,eimag;

main()
{double a,b,c,d,x,y,m,k,e,omp,eh,eta,g2,g3;
struct complex z,p,pp,etap,ans;
struct complex v,qq,ans3,ans1,ans2;int i,j;double zz;
struct complex ee,f,g;
double r,s,t,sig=1.e-4;
double ts,tc,td,tn,u,alpha;
double eps,n,phi;

infinite_loop
	{
	printf(" E(u) enter mag u, theta(deg), m(==0. to quit)\n");
	scanf("%le%le%le",&a,&b,&d);
	c= b*pi/180.;
	if(d==0.)break;
	else printf(" for m=%le theta(rad)=%le,|u|=%le\n",d,c,a);
	CMPLX(z,a*cos(c),a*sin(c));
	ce(&z,d,&p);printc(&p);
	}
infinite_loop
	{
	printf(" K,E(m) enter m(==0. to quit)\n");
	scanf("%le",&a);
	if(a==0.)break;
	tek(0,a,&c,&d);
	printf(" K,E= %le %le\n",c,d);
	if(kimag!=0.)printf(" imag parts %le %le\n",kimag,eimag);
	}

while(1)
	{
	a=0.;b=0.;
	printf(" for tef(incomplete): enter m (-1 to quit) ");
	scanf("%le",&x);m=x; if(m==-1.)break;
	printf("enter phi");scanf("%le",&alpha);
	/*if(alpha==0. && m==0.)break;*/
	tef(alpha,m,1.e-5,&a,&b);
	printf(" at phi=%le,m=%le,F=%le E=%le Z=%le\n"
		,alpha,m,a,b,jzeta(alpha,m));
	}

while(1)
	{
	/*alpha= (j+1)*pi*15./180.; m=sin(alpha);m*=m;*/
	a=0.;b=0.;
	printf(" tef: enter m -1. to quit ");
	scanf("%le",&x);m=x;
	if( m==-1.)break;
	tef(pi/6.,m,1.e-3,&a,&b);printf(" at pi/6=30,m=%g,%le %le\n",m,a,b);
	}
printf(" complete integrals\n");
for(j=0;j<16;j++)
	{
	m= (j+1)*.05;
	tek(0,m,&a,&b);printf(" %le %le %le\n",m,a,b);
	}
printf(" incomplete, 30deg \n");
for(j=0;j<4;j++)
	{
	alpha= (j+1)*pi*15./180.; m=sin(alpha);m*=m;a=0.;b=0.;
	printf(" alpha=%f %g ",alpha,m);
	tef(pi/6.,m,1.e-3,&a,&b);printf(" %le %le\n",a,b);
	}
printf(" jzeta\n");
printf(" %le \n",jzeta( pi/6.,.25));
printf(" jacobian elliptic\n");
jef(.2,.81,&a,&b,&alpha);
printf(" .2|.81 sn=%le cn=%le dn=%le %le\n",a,b,alpha,jtheta);
jef(.2,.19,&a,&b,&alpha);
printf(" .2|.19 sn=%le cn=%le dn=%le jtheta=%le\n",a,b,alpha,jtheta);
jef(.75342,.7,&a,&b,&alpha);
printf(" .75342|.7 sn=%le cn=%le dn=%le %le\n",a,b,alpha,jtheta);
jef(.6,.36,&a,&b,&alpha);
printf(" .6|.36 sn=%le cn=%le dn=%le %le\n",a,b,alpha,jtheta);
alpha=heuman(pi/6.,.25);
printf(" heuman=%le\n",alpha  );

eps=.3333333;/* 30 degrees */
u= eps*1.68575;
neville(u,.25,&ts,&tc,&td,&tn);
printf(" neville s,c,d,n %le %le %le %le\n",ts,tc,td,tn);



while(1){
printf(" e3: enter n phi m  all 0 to end.\n");
scanf("%le%le%le",&a,&b,&c);
n=a;phi=b;m=c;
if(a==0. && b==0. && c==0.)break;
printf(" %g %f %f %le\n",n,phi,m,e3(n,phi,m,1.e-4));
}




while(1)
	{
	printf(" Weierstrass:enter z=x+iy,g2,g3 quit z==0");
	scanf("%le%le%le%le",&a,&b,&c,&d);
	if(a==0.)break; z.x=a;z.y=b; x=c;y=d; g2=c;g3=d;
	weier(&z,x,y,&p,&pp,&m,&k,&e,&omp,&eta,&eh);
	printf(" P and P':");printc(&p);printf(" ");printc(&pp);printf("\n");
	printf(" m= %le K=%le E=%le omega(2) %le e or h2 %le\n",m,k,e,omp,eh);
	sigma(&z,&ans,k,omp,eta,m,g2,g3);
	printf(" sigma= %le %le\n",ans.x,ans.y);
	zetaw(&z,&ans,k,omp,eta,m,g2,g3);
	printf("zetaw %le %le\n",ans.x,ans.y);
	}
printf(" m given q\n");
for(i=1;i<21;i++)
	{
	zz= i/20.;qq.y=0.;qq.x=zz;emf(&qq,&ans);
	printf(" q=%le m=%le or %le %le\n",zz,mq(zz),ans.x,ans.y);
	}
for(i=1;i<21;i++)
	{
	zz= i/100.;
	printf(" q=%le m=%le %le\n",zz,mq(zz),mq(zz));
	}

printf(" q for .9997 %le\n",q(.9997));

while(1){
printf(" theta:enter v,q both complex");scanf("%le%le%le%le",&a,&b,&c,&d);
if(a==0. && b==0. && c==0. && d==0.)break;
CMPLX(v,a,b); CMPLX(qq,c,d);ctheta(&v,&qq,&ans1,&ans2,&ans3,&ans);
printf(" ans=%le %le %le %le %le %le %le %le\n"
,ans.x,ans.y,ans3.x,ans3.y,ans1.x,ans1.y,ans2.x,ans2.y);
}

quartic(0.,0.,0.,-1.,&z,&ee,&f,&g);
printf(" quartic %le %le %le %le %le %le %le %le\n"
	,z.x,z.y,ee.x,ee.y,f.x,f.y,g.x,g.y);
while(1)
{
printf(" enter real m, complex z for cef");scanf("%le%le%le",&a,&b,&c);
if(a==0.&& b==0. &&c==0.)break;
m=a;z.x=b;z.y=c; cef(m,&z,&ee,&f,sig);
printf(" e %le %le f %le %le\n",ee.x,ee.y,f.x,f.y);
}
while(1)
{printf(" enter complex z, real g2 g3 for invw test\n");
scanf("%le%le%le%le",&a,&b,&c,&d);
g2=c;g3=d; z.x=a;z.y=b;
invp(g2,g3,&z,&f);
printf(" answer %le %le\n",f.x,f.y);
}
}