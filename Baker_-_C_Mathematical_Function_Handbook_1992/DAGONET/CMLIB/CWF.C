/*  Coulomb Wave Functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define min(a,b) ((a)<(b)?(a):(b))
#define tol 1.e-6
#define MAXIT 500

static double g0;

double coulombf(eta,rho,l) int l;double eta,rho;
{
/* flp is derivative */
int i,small;
double exp(),pow(),log(),term,sum,sqrt(),c,cf,c0,power,a,aold,aolder,gl,fl;
double theta,psi,d,qp,sump,flp,termp,rhol;
double sa,saold,saolder,pl,rl,terms;
if(rho>35.)
	{cwfa(l,eta,rho,&fl,&gl,&flp,&qp);
	if(!l)g0=gl;
	return fl;
	}
term=2.*pi*eta;
if( eta==0.)
	{
	/*if(!l)return sin(rho);*/
	return (rho)*sjn(rho, l);
	}
if(eta<1.e-10)c0=1.;
else c0=sqrt(term/(exp(term)-1.));
c=c0;
for(i=1;i<=l;i++)
     c*=sqrt( eta*eta+i*i)/(i*(2*i+1));
/*printf(" c0=%e\n",c0);*/
power=rho; aolder=1.;aold=eta/(1.+l);sum=aolder+aold*power;
sump= (l+1.+ aold*(l+2.)*power);
saolder=1.;saold=0.;
/*d= 1./(c*(2*l+1.));*/
pl= c*c/(c0*c0)*eta*(2*(2*l+1));
i=2;
if(!l)
	{sa= (2.*eta*saold-saolder-pl*(2*i-1)*aold)/((i-l-1)*(i+l));
	psi=1.+ saold*rho+sa*rho*rho;
	saolder=saold;saold=sa;
	}
for(small=0,i=3;i<MAXIT;i++)
	{
	power*=rho;
	a= (2.*eta*aold-aolder)/((double)((l+l+i)*(i-1)));
	if(!l)
		{
		sa= (2.*eta*saold-saolder-pl*(2*i-1)*a)/((i-l-1)*(i+l));
		terms=sa*power*rho;
		psi+=terms;
		}
	term=power*a;
	sum+=term;
	termp=term*(i+l);
	sump+=termp;
	/*printf(" sum, term, %le %le %d\n",sum,term,i);*/
	if( abs(sum)*tol>abs(term)&&abs(sump)*tol>abs(term))
		{if(!small)small=1;
		if(small)break;/* two in a row*/
		}
	else small=0;
	aolder=aold;aold=a;
	saolder=saold;saold=sa;
	}
if(l)	rhol=pow(rho,(double)(l));
else	rhol=1.;
fl= c*sum*rho*rhol;
flp= c*sump*rhol;
if(rho>34.5)
	{
	cwfa(l,eta,rho,&saold,&gl,&saolder,&qp);
	theta= ((rho-34.5)/(35.-34.5));
/*printf(" theta=%e %e %e\n",theta,fl,saold);*/
	fl= fl*(1.-theta)+saold*(theta);
	flp=flp*(1.-theta)+saolder*(theta);
	if(!l)g0=gl;
	}
/*printf(" flp=%e\n",flp);*/
return fl;
}

static double etaparm,rhoparm;
static int lparm;

double gl(q) double q;
{double x,y,z,pow(),exp(),tanh(),atan(),sin(),tan(),t,cos(),jacob,ans;
if(q<.5*pi)t=tan(q);
else t=1.e12;
x=tanh(t);y=1.-x*x;z=1.;
if(lparm>0){z=1.+t*t;if(lparm>1)z=pow(z,(double)lparm);}
jacob=cos(q);
jacob=1./(jacob*jacob+1.e-8);
if(lparm>0)y=pow(y,1.+lparm);
ans= (z*exp(-rhoparm*t+2.*etaparm*q)
	-y*sin(rhoparm*x-2.*etaparm*t))*jacob;
return ans;
}

double coulombg(eta,rho,l) int l;double eta,rho;
{
int i; double term,coef,c0,c,exp(),pow(),gamma(),adsimp(),as,sqrt(),btm,top,delta;
etaparm=eta;rhoparm=rho,lparm=l;
/*if(!l && eta<1.)
	{coef=coulombf(eta,rho,l);return g0;
	}
*/
term=2.*pi*eta;
if(eta<1.e-10)c0=1.;
else c0=sqrt(term/(exp(term)-1.));
c=c0;
for(i=1;i<=l;i++)
     c*=sqrt( eta*eta+i*i)/(i*(2*i+1));
coef=exp(-pi*eta)*pow(rho,1.+l)/(gamma(2*l+2.)*c);
as=0.;
btm=0.;top=.025;delta=.05;
for( i=0;btm<.5*pi;i++)
	{as+=adsimp(btm ,top ,tol,gl);btm=top;
	top= min(.5*pi,top+delta);}
/*as=adsimp(0.,.02,tol,gl)+adsimp(.02,.1,tol,gl)+adsimp(.1,.25*pi,tol,gl)
+adsimp(.25*pi,.5*pi,tol,gl);*/
return as*coef;
}

cwfa(l,eta,rho,fl,gl,flp,glp) int l; double eta,rho,*fl,*gl,*flp,*glp;
{
struct complex carg,cans,cdummy;
double c,theta,sigma,f,g,fs,gs,fold,gold,d,ct,st,log(),argmt(),
	fsold,gsold,a,b,fsum,gsum,fssum,gssum,cos(),sin();
int k;
fold=1.;gold=fsold=0.;gsold=1.-eta/rho;
fsum=fold;gsum=gold;gssum=gsold;fssum=fsold;d=eta*eta+l*(l+1);
for( k=1;k<100;k++)
	{
	c=1./((2*k+2)*rho);
	a=(2*k+1)*eta*c;
	b=(d-k*(k+1))*c;
	f=a*fold-b*gold;
	g=a*gold+b*fold;
	fs=a*fsold-b*gsold-f/rho;
	gs=a*gsold+b*fsold-g/rho;
	if(k>1)
		{
		if( abs(fsum)*tol>abs(f) || abs(f)>abs(fold))break;
		if( abs(gsum)*tol>abs(g) || abs(g)>abs(gold))break;
		}
	fsum+=f;gsum+=g;fssum+=fs;gssum+=gsum;
	gold=g;fold=f;fsold=fs;gsold=gs;
	}
carg.x=1.+l;carg.y=eta;cgamma(&carg,&cans,&cdummy);
sigma=argmt(cans.y,cans.x);
theta= rho-eta*log(2.*rho)-.5*pi*l+sigma;
ct=cos(theta);st=sin(theta);
*fl=gsum*ct+fsum*st;
*gl=fsum*ct-gsum*st;
*flp=gssum*ct+fssum*st;
*glp=fssum*ct-gssum*st;
return 0;
}
