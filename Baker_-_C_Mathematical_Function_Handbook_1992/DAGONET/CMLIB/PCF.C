/*  Parabolic cylinder functions 
and related confluent hypergeometric functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

upcf	U
vpcf	V
wpcf	W: Associated pcf (for x, -x independent solutions)

uses confluent hypergeometric function

*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"



double wpcf(a,xx) double a,xx;
{
int i,twoi,itmax=100;
struct complex ag,loggam,gam;
double g1,g3,y1,y3,term1,term3,x,tol=1.e-7,wa(),sqrt();
double z,oold,oolder,eold,eolder,new,f1,f3,sign,wairy();
x=abs(xx);
if(a>4.) return wairy(a,xx);
if( x>8. && x> abs(a)*4.) return wa(a,xx);

if(xx>=0.)
	{
	sign=-1.;
	}
else
	{
	sign=1.;
	}
ag.x=.25;ag.y=.5*a;
cgamma(&ag,&gam,&loggam);
g1= cabs(gam);
ag.x=.75;ag.y=.5*a;
cgamma(&ag,&gam,&loggam);
g3= cabs(gam);
z=x*x;
oolder=eolder=1.;
eold=oold=a;
f1=.5*z;
f3=.1666666*z;
y1=1.+a*f1;
y3=1.+a*f3;
for (i=2;i<itmax;i++)
	{
	twoi=i<<1;
	new= a*eold-.25*(twoi-2)*(twoi-3)*eolder;
	f1*= z/((twoi)*(twoi-1));
	term1=f1*new;
	y1+=term1;
	eolder=eold;
	eold=new;
	new= a*oold-.25*(twoi-2)*(twoi-1)*oolder;
	f3*= z/((twoi+1)*(twoi));
	term3=f3*new;
	y3+= term3;
	oolder=oold;
	oold=new;
	if( abs(term1)<tol && abs(term3)<tol)break;
	}
y3*=x;
f1=sqrt(g1/g3);
f3=sign*sqrt(2.*g3/g1);
return  sqrt(sqrt(2.))*.5*(f1*y1+f3*y3);
}

double wa(a,xx) double a,xx;
{
int i,itmax=6;
struct complex ag,loggam,ans,ans1;
double arg,sqrt(),log(),cos(),sin(),x,y,phi,s1,s2,k,exp(),
term1,term2,factor,tol=1.e-6,old1,old2;
printf(" entered wa\n");
if(a<=0.) return errorcode;
ag.x=.5; ag.y=a;
x=abs(xx);
y=exp(pi*a);
k= sqrt(1.+y*y)-y;
cgamma(&ag,&ans,&loggam);
phi=loggam.y;
arg= .25*x*x-a*log(x)+.25*pi+.5*phi;
y=.5/(x*x);
factor=y;
s1=1.;
s2=0.;
old1=old2=1.e10;
for(i=1;i<itmax;i++)
	{
	ag.x=.5+((double) (i<<1));ag.y=a;
	cgamma(&ag,&ans1,&loggam);
	CDIV(loggam,ans1,ans);
	switch (i%4)
		{
		case 1:
			term1=loggam.y;
			term2=-loggam.x;
			break;
		case 2:
			term1=-loggam.x;
			term2=-loggam.y;
			break;
		case 3:
			term1=-loggam.y;
			term2=loggam.x;
		default:
			term1=loggam.x;
			term2=loggam.y;
		}
	term1*=factor;
	term2*=factor;
	if(abs(term1)>abs(old1) || abs(term2)>abs(old2)) break;
	s1+=term1;
	s2+=term2;
	if( abs(term1)<tol && abs(term2)<tol)break;
	if( abs(term1)<tol*abs(s1) && abs(term2)<tol*abs(s2))break;
	factor*=y/(i+1);
	old1=term1;
	old2=term2;
	}
if(xx>=0.)
	{return sqrt(2.*k/(x))*(s1*cos(arg)-s2*sin(arg));}
/*else*/
	return sqrt(2./(k*x))*(s1*sin(arg)+s2*cos(arg));
}

#define acosh arccosh

double wairy(a,xx) double a,xx;
{
double sqrt(),pow(),acos(),acosh(),bi(),ai(),xi,t,tau,th,x,exp(),airy,ex,c;
/* a>>0 */
xi= abs(xx/(2.*sqrt(a)) ) ;
if(xi<1.)
	{
	th= .25*(acos(xi)-xi*sqrt(1.-xi*xi));
	t=-pow(1.5*th,.6666666);
	}
else
	{
	th= .25*(xi*sqrt(xi*xi-1.)-acosh(xi));
	t=pow(1.5*th,.6666666);
	}
t*=pow(4.*a,.666666);
if(t!=0.)c=sqrt(pi*sqrt(t/((4.*a)*(xi*xi-1.))));
else c=sqrt(pi/sqrt(4.*a));
ex=exp(-.5*pi*a);
if(xx>0.)
	return c*ex*bi(-t);
return 2.*c/ex*ai(-t);

}

double upcf(nu,x)double nu,x;
{
double u(),m(),mult=4.,sqrt(),gamma();
double a,b,c,d;

if( x> 4. && x> abs(nu)*mult)
	{
	a=.5/(x*x);
	d= exp(-.25*x*x)/pow(x,nu+.5)*
		(1.+(nu+.5)*(nu+1.5)*a*((nu+2.5)*(nu+3.5)*a*.5-1.));
	printf(" asymptotic u=%e\n",d);
	}

if(x>0.)
return exp(-.25*x*x)/pow(2., .25+.5*nu)*u(.25+.5*nu,.5,x*x*.5);
/*else*/
a= sqrt(pi)*exp(-.25*x*x)/pow(2.,.25+nu*.5);
b=gamma(.75+.5*nu);
c=gamma(.25+.5*nu);
d= a*
	(m(.5*nu+.25,.5,.5*x*x)/b
	-x*sqrt(2.)*m(.5*nu+.75,1.5,.5*x*x)/c);
return d;
}

double vpcf(nu,x)double nu,x;
{
double a,b,c,d,mult=4.;
if( x> 4. && x> abs(nu)*mult)
	{
	a=.5/(x*x);
	d= exp(.25*x*x)*pow(x,nu-.5)*sqrt(2./pi)*
		(1.+(nu-.5)*(nu-1.5)*a*((nu-2.5)*(nu-3.5)*a*.5+1.));
printf(" asymptotic v=%e\n",d);
	}
a= gamma(nu+.5)/pi;
b=sin(pi*nu);
c=upcf(nu,x);
d=upcf(nu,-x);
return a*(b*c+d);
}

/* whittaker functions M and W */

double mwhit(k,mu,z) double k,mu,z;
{
double m();
return m(.5+mu-k,1.+2.*mu,z)*exp(-z*.5)*pow(z,.5+mu);
}

double wwhit(k,m,z) double k,m,z;
{
double u();
return u(.5+m-k,1.+2.*m,z)*exp(-z*.5)*pow(z,.5+m);
}

/* simple versions of confluent hypergeometric functions M and U
for real arguments*/

double m(a,b,x)
double a,b,x;
{
double sum,term,newterm,count,tol=1.e-6,pow(),exp(),gamma(),oldterm;
int top=40,i,k,itmax=40;

if(x==0.)return 1.;
if((b-(double)((int)b))==0.  && b<0. &&
	(a-(double)((int)a))==0. && b<a)return(errorcode);
if(x> 10.)
	{
	sum=1.;term=1.;oldterm=1.e30;
	for (i=0;i<itmax;i++)
		{
		k=i+1;
		term*=(b-a+i)*(k-a)/(x*(k));
		newterm=abs(term);
		if(newterm>oldterm)break;
		oldterm=newterm;
		sum+=term;
		if((abs(term)<tol||abs(term)<tol*abs(sum)))break;
		}
	return exp(x)*pow(x,a-b)*gamma(b)/gamma(a)*sum;
	}
for(i=0,term=sum=1.;i<top;i++)
	{
	count=(double)(i+1);
	term *=(x*(a+i)/((b+i)*count));
	if( abs(term)< tol*abs(sum))return sum;
	sum+=term;
	}
return(sum);
}

double u(a,b,x)
double a,b,x;
{
double pow(),sin(),fac(),m(),eps=1.e-4,gamma(),tol=1.e-6,exp(),cos();
double sum,term,count,ret,factor;
int i,top=40;

/* wont work for integer b*/
 if(abs(x)<5.)
	{
	help:
/*printf(" after help x, a b %f %f %f\n",x,a,b);*/
if((b-(double)((int)b))==0. ) {b+=eps;/* avoid sin(pi*b)=0*/}
	if(x>0.)
	return pi/(sin(pi*b))*(m(a,b,x)/(gamma(1.+a-b)*gamma(b))
	-m(a+1.-b,2.-b,x)/(gamma(a)*fac(1.-b))*pow(abs(x),1.-b));
	sum=-(x);
	/* factor= exp i pi (1-b) */
	factor= cos( pi*(1.-b));/*neglect imaginary part sin(), if any for now*/
	return pi/(sin(pi*b))*exp(-sum)*(m(b-a,b,sum)/(gamma(1.+a-b)*gamma(b))
	-m(1.-a,2.-b,x)*factor/(gamma(a)*gamma(2.-b))*pow(sum,1.-b));
	}
/*else*/
/*printf(" asympt u a=%f b=%f x=%f\n",a,b,x);*/
for (i=0,term=1.,sum=1.;i<top;i++)
	{
	count=(double)(i+1);
	factor=(a+i)*(a-b+count)/(-x*count);
	if(abs(factor)>1.)
	/* asympt. series invalid once terms increase*/
		{
	/*	printf(" term=%d sum=%f term=%f,factor=%f\n",i,sum,term,factor);*/
		if (abs(term/sum)<1.e-3)return sum/pow(x,a);
/*		printf(" getting help,would be %f\n", sum/pow(x,a));*/
		goto help;
		}
	term*=factor;
	sum +=term;/*sum must be nonzero here */
	if(abs(term/sum)<tol)return (sum/pow(x,a));
	}
/*printf(" asympt non converge %f %f %f\n",factor,term,sum);*/
return sum/pow(x,a);/*keep desmet happy*/
}
