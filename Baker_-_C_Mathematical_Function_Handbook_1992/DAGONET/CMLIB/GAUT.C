/*  Legendre functions for |x|>1
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
based upon CALGO 259 by Walter Gautschi
		
leg1(x,a,nmax,p) double x,p[];int a,nmax;  Pa,n(x) n=0 to nmax
leg2(x,m,nmax,d,q) double x,q[];int m,nmax,d; Qm,n(x)  n=0 to nmax
leg3( x,n,mmax,d,q) int n,mmax,d; double x,q[]; Qm,n(x) m=0 to mmax
legend1(x,alpha,nmax,d,p)      double x,alpha,p[]; int d,nmax; Pa,n=0 nmax
legend2(x,a,m,nmax,d,p) double p[],a,x;int m,nmax,d;P a+n,m n=0 to nmax
conical(x,tau,nmax,d,p) double x,tau,p[];int nmax,d;
toroidal(x,m,nmax,d,q)   int m,d,nmax; double x,q[];

not based upon CALGO 259:
arccosh(x)	arc hyperbolic cosine
double conicalt(theta,tau) double theta,tau;conical for |x|<1
double conicala(x,tau) double x,tau; asymptotic conical

*/
#include <stdio.h>
#include <alloc.h>
#include "cmlib.h"
#include "protom.h"
/*
   #define malloc farmalloc
   #define free   farfree
*/

leg1(x,a,nmax,p) double x,p[];int a,nmax;
{/* Pa,n(x) n=0 to nmax x>1*/
double *rr, x1,sum,r,s,pow(),sqrt();       int n,nstop;
if( x<1. || a<0 || nmax<0)
	{fprintf(stderr," leg1 bad argument\n");return ierrorcode;}
if(x==1. ||a==0){p[0]=1.;for(n=1;n<=nmax;n++)p[n]=0.;return 0;}
rr=malloc( (nmax+1)* sizeof(double) );
if(rr==NULL)
	{fprintf(stderr," leg1 memory allocation error\n");return ierrorcode;}
for(n=a+1;n<=nmax;n++) p[n]=0.;
x1=sqrt(x*x-1.);
sum= pow( (x+x1),(double)a);x1=2.*x/x1;
r=s=0.;
for(n=a;n>=1;n--)
	{
	r= (a+1-n)/(n*x1+r*(n+a+1));
	s=r*(2.+s);
	if(n<=nmax) rr[n-1]=r;
	}
p[0]=sum/(1.+s);
nstop= (nmax<=a)? nmax-1:a-1;
for(n=0; n<=nstop;n++)p[n+1]=(n+a+1.)*rr[n]*p[n];
free(rr);
return 0;
}


leg2(x,m,nmax,d,q) double x,q[];int m,nmax,d;
{/* qn,m(x) n=0 to nmax x>1*/
int n,nu,p,flag;
double sqrt(),pow(),log(),x1,q0,q1,q2,epsilon,r,*rr,*qapprox;
if( x<1. || m<0 || nmax<0)
	{fprintf(stderr," leg2 bad argument\n");return ierrorcode;}
rr=malloc( (nmax+1)* sizeof(double) );
qapprox=malloc( nmax* sizeof(double) );
if(rr==NULL||qapprox==NULL)
	{fprintf(stderr," leg2 memory allocation error\n");return ierrorcode;}
x1=sqrt(x*x-1.);
q1= .5*log((x+1.)/(x-1.));
if(!m)	q[0]=q1;
else
	{
	q2=-1/x1;x1=2.*x/x1;
	for(n=1;n<m;n++)
		{
		q0=q1;q1=q2;
		q2=-n*x1*q1-n*(n-1)*q0;
		}
	q[0]=q2;
	}
for(n=0;n<=nmax;n++)qapprox[n]=0.;
epsilon=.5*pow(10.,(double)-d);
nu=20+ 1.25*nmax;
while(1)
	{r=0.;
	flag=1;
	for(n=nu; n>=1;n--)
		{
		r=(n+m)/(((n<<1)+1)*x-(n-m+1.)*r);
		if(n<=nmax)rr[n-1]=r;
		}
	for(n=0;n<nmax;n++)q[n+1]=rr[n]*q[n];
	for(n=0;n<=nmax;n++)
		{
		if( abs(q[n]-qapprox[n])>epsilon*abs(q[n]) ){flag=0;break;}
		}
	if(flag)break;
	for(p=0;p<=nmax;p++)qapprox[p]=q[p];
	nu+=10;
	if(nu> 400){fprintf(stderr," leg2 cannot meet tolerence req.\n");break;}
	}
free(rr);free(qapprox);return 0;
}


leg3( x,n,mmax,d,q) int n,mmax,d; double x,q[];
{/* qn,m m=0 to mmax*/
int m; double x1,*q1,sqrt();
if(n<0 ||x<1. ||mmax<0)
	{fprintf(stderr," leg3 argument error\n");return ierrorcode;}
q1=malloc((n+1)* sizeof(double) );
if(q1==NULL)
	{fprintf(stderr," leg3 memory allocation error\n");return ierrorcode;}
leg2(x,0,n,d,q1); q[0]=q1[n];
x1=2.*x/sqrt(x*x-1.);
if(mmax>0)
	{leg2(x,1,n,d,q1); q[1]=q1[n];}
for(m=1;m<mmax;m++)q[m+1]=-m*x1*q[m]-(m+n)*(m-n-1.)*q[m-1];
free(q1);
return 0;
}


legend1(x,alpha,nmax,d,p)      double x,alpha,p[]; int d,nmax;
{/* p alpha n n=0 nmax x>1*/
int n,nu,m,flag; double a,epsilon,x1,sum,c,r,s,*papprox,*rr,pow(),sqrt();
if(x<1. ||nmax<0){fprintf(stderr," legend1 argument error\n");return 0;}
if(x==1.)
	{
	p[0]=1.;for(n=1;n<=nmax;n++)p[n]=0.;return 0;
	}
rr= malloc ( (nmax+1)* sizeof(double) );
papprox= malloc ( (nmax+1)* sizeof(double) );
if(rr==NULL||papprox==NULL)
	{fprintf(stderr," legend1 memory allocation error\n");return ierrorcode;}
a=(alpha<.5)? -alpha-1.: alpha;
for(n=0;n<=nmax;n++)papprox[n]=0.;
epsilon=pow(10.,(double) -d)*.5;
x1=sqrt(x*x-1.); sum=pow(x+x1, a);x1=2.*x/x1;
nu=20+ (37.26+.1283*(a+38.26)*x)*nmax/(37.26+.1283*(a+1.)*x);
while(1)
	{
	flag=1;
	r=s=0.;
	for(n=nu;n>=1;n--)
		{
		r=(a+1.-n)/(n*x1+(n+a+1.)*r);s=r*(2.+s);
		if(n<=nmax)rr[n-1]=r;
		}
	p[0]=sum/(1.+s);
	for(n=0;n<nmax;n++)p[n+1]=rr[n]*p[n];
	for(n=0;n<=nmax;n++)
		{
		if(abs(p[n]-papprox[n])>epsilon*abs(p[n])){flag=0;break;}

		}
	if(flag)break;
	for(m=0;m<=nmax;m++)papprox[m]=p[m];nu+=10;
	if(nu> 200){fprintf(stderr," legend1 cannot meet tolerence req.\n");break;}
	}
c=1.;
for(n=1;n<=nmax;n++)
	{c*=(a+n);p[n]*=c;
	}
free(rr);free(papprox);
return 0;
}

legend2(x,a,m,nmax,d,p) double p[],a,x;int m,nmax,d;
{/*p a+n,m n=0 nmax x>1*/
double *p1,denom; int n;
if(m<0 ||x<1. ||nmax<0){fprintf(stderr," legend2 argument error\n");return ierrorcode;}
p1= malloc ((nmax+1) * sizeof(double) );
if(p1==NULL){fprintf(stderr," legend2 memory allocation error\n");return ierrorcode;}
legend1(x,a,m,d,p1); p[0]=p1[m];
if(nmax>0)
	{
	legend1(x,a+1.,m,d,p1);
	p[1]=p1[m];
	}

for(n=1;n<nmax;n++)
	{
	denom=(n-m+1)+a;
	if(denom==0.)
			{fprintf(stderr,"legend2 zero denominator n=%d\n",n);
			p[n+1]=errorcode;
			return 1;}
	p[n+1]=(((n<<1)+2.*a+1.)*x*p[n]-(n+a+m)*p[n-1])/(denom);
	}
return 0;
}

#define numax 1000

conical(x,tau,nmax,d,p) double x,tau,p[];int nmax,d;
{/* p -.5+i*tau,n (x) n=0 to nmax x>1 */
int n,nu,m,flag;
double epsilon,t,x1,x2,sum,lambda1,lambda2,lambda,r,s,
	*rr,*papprox,pow(),sqrt(),log(),cos(),aux1,aux2,aux3;
if(x<1. || nmax<0){ fprintf(stderr," conical argument error\n");return 0;}
if(x==1.)
	{p[0]=1.;for(n=1;n<=nmax;n++)p[n]=0.;return 0;}
rr=malloc( (nmax+1)* sizeof(double) );
papprox=malloc( (nmax+1)* sizeof(double) );
if(rr==NULL||papprox==NULL)
	{fprintf(stderr," conical memory allocation error\n");return 0;}
t=tau*tau;
for(n=0;n<=nmax;n++)papprox[n]=0.;
epsilon=.5*pow(10.,(double)-d);
x1=sqrt(x*x-1.);x2=x+x1;
sum= cos(tau*log(x2))/sqrt(x2); x1=2.*x/x1;
nu=30+ ((1.+(.14+.0246*tau)*(x-1.))*nmax);
while(1)
	{
	flag=1;n=2;
	lambda1=1./(.25+t);
	lambda2=(3.-4.*t)/((.25+t)*(2.25+t));
	L1:aux2=tau/n; aux1=1.+.5/n;
	lambda= (1.+1./n)*(2.*lambda2-lambda1)/(aux1*aux1+aux2*aux2);
	if(n<nu)
		{lambda1=lambda2;lambda2=lambda;n++;goto L1;}
	r=s=0.;
	L2:
	aux1= 1.-.5/n;aux2=tau/n;
	aux3=(aux1*aux1+aux2*aux2);
	r=-aux3/(x1+(1.+1./n)*r);
	s=r*(lambda2+s);
	if(n<=nmax) rr[n-1]=r;
	lambda1=lambda2;
	aux1=1.+.5/n;
	lambda2=2.*lambda2-(aux2*aux2+aux1*aux1)*lambda/(1.+1./n);
	lambda=lambda1;n--;
	if(n>=1)goto L2;
	p[0]=sum/(1.+s);
	for(n=0;n<nmax;n++)p[n+1]=rr[n]*p[n];
	for(n=0;n<=nmax;n++)
		{
		if(abs(p[n]-papprox[n])>epsilon*abs(p[n])){flag=0;break;}
		}
	if(flag)break;
	for(m=0;m<=nmax;m++)papprox[m]=p[m];nu+=60;
	if(nu>numax){fprintf(stderr," conical did not meet tolerence req.\n");break;}
	}

t=1.;
for(n=1;n<=nmax;n++) {t*=n;p[n]*=t;}
free(rr);free(papprox); return 0;
}


toroidal(x,m,nmax,d,q)   int m,d,nmax; double x,q[];
{/* q -.5+n,m x n=0 nmax x>1*/
double *rr,*qapprox,pow(),sqrt(),epsilon,r,s,sum,c,x1;int flag,n,nu;
if(x<=1. || nmax<0){ fprintf(stderr," toroidal argument error\n");return ierrorcode;}
rr=malloc( (nmax+1)* sizeof(double) );
qapprox=malloc( (nmax+1)* sizeof(double) );
if(rr==NULL||qapprox==NULL)
	{fprintf(stderr," toroidal memory allocation error\n");return ierrorcode;}
for(n=0;n<=nmax;n++)qapprox[n]=0.;
epsilon= pow(10.,(double)-d)*.5;
c=2.2214414691;
if(m>=0)
	for(n=0;n<m;n++)c*= -(n+.5);
else
	for(n=0;n>=m+1;n--)c/=(.5-n);
sum=c/sqrt(x-1.);if(m)sum*=pow((x+1.)/(x-1.),.5*m);x1=2.*x;
nu=20+ (1.15+(.0146+.00122*m)/(x-1.))*nmax;
while(1)
	{flag=1;
	r=s=0.;
	for(n=nu;n>=1;n--)
		{
		r=(n+m-.5)/(n*x1-r*(n-m+.5));
		s=r*(2.+s);
		if(n<=nmax)rr[n-1]=r;
		}
	q[0]=sum/(1.+s);
	for(n=0;n<nmax;n++)q[n+1]=q[n]*rr[n];
	for(n=0;n<=nmax;n++)
		{
		if(abs(q[n]-qapprox[n])>epsilon*abs(q[n])){flag=0;break;}
		}
	if(flag)break;
	for(n=0;n<=nmax;n++)qapprox[n]=q[n];nu+=10;
	if(nu>numax){fprintf(stderr," toroidal did not meet tolerence req.\n");break;}
	}
free(rr);free(qapprox);return 0;
}

#define itmax 60
#define tol 1.e-9

double conicalt(theta,tau) double theta,tau;
{/* conical for -1<cos(theta)<=1 as argument. series*/
double s,sum,t,term,sin(); int j,i,k;
if(theta<0. || theta>=pi)return errorcode;
s=sin(.5*theta);s*=s;sum=1.;term=1.;t=tau*tau;
for(i=1;i<itmax;i++)
	{
	j=(i<<1);k=j-1;
	term*=(t+.25*k*k)*s/(i*i);
	sum+=term;
	if(abs(term)<tol*abs(sum))return sum;
	}
fprintf(stderr," no convergence conical x<1\n");
return sum;
}

double arccosh(x) double x;
{double log(),sqrt();
return  log( x+sqrt(x*x-1.));
}

double conicala(x,tau) double x,tau;
{double b,acos(),sqrt(),exp(),sin(),cos(),sinh();
/* P .5-i*tau (x) for large tau */
if(abs(x)>1.)
	{
	b=arccosh(x);
	return sqrt(2./(pi*tau*sinh(b)))*sin(tau*b+.25*pi);
	}
b=acos(x);
return exp(tau*b)/sqrt(2.*pi*tau*sin(b));
}
