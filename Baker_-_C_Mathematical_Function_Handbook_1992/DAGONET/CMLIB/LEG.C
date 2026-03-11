/* Legendre and Associated Legendre functions
Copyright 1991 by Louis Baker. All rights reserved.

	real arguments:
double plm(l,m,x) int l,m;double x;
double pl0(l,x) double x;
double ql0(l,x) double x;
double qlm(l,m,x) double x; int l,m;
double pmunu(mu,nu,x) double nu,mu,x;
double qmunu(mu,nu,x) double nu,mu,x;
double pnu(x,nu)double x,nu;
	imaginary arguments:
double pli(l,z) int l; double z;
double qli(l,z) int l; double z;
	real=1 for real arguments, 0 imaginary:
double legendrea(m,n,x,real) int m,n,real;double x;
	return a table in q[] (ratios in r[]) for given m, n=0,nmax
qleg(m,nmax,x,real,r,q) int m,nmax,real; double x,r[],q[];
double qnu(x,nu,real)int real;double x,nu;

legendrea based upon CALGO 47,qleg CALGO 62 by J. R. Herndon
*/


#include "cmlib.h"
#include "complex.h"
#include "protom.h"
#include <stdio.h>

double plm(l,m,x) int l,m;double x;
{/* m<=l*/
int i;
double pmm,coef,factor,sqrt(),pl0();
/*if(abs(x)>1.)return errorcode;*/
if(m==0 && abs(x)>1.)return pl0(l,x); /* do not use pl0 off cut*/
if( abs(x)>1.)
	{
	return ((m-1-l)*x*plm(l,m-1,x)-(m+l-1)*plm(l-1,m-1,x))/sqrt(x*x-1.);
	/* will ultimately get to pl0 call*/
	}
pmm=1.;/* |x|<=1 */
if(m>0)
	{
	coef=sqrt((1.-x)*(1.+x));
	factor=1.;
	for(i=1;i<=m;i++)
		{
		pmm*=-factor*coef;
		if(pmm==0.)break;/* x=1*/
		factor+=2.;
		}
	}
if(l==m)return pmm;
coef= x*( (m<<1)+1.)*pmm;
if(l==(m+1))return coef;
for(i=m+2;i<=l;i++)
	{factor= (x*((i<<1)-1)*coef-(i+m-1)*pmm)/(i-m);
	pmm=coef;
	coef=factor;
	}
return factor;
}

double pl0(l,x) int l; double x;
{int i;
double old,older,new;
if(l==0)return 1.;
if(l==1)return x;
old=x;older=1.;
for(i=2;i<=l;i++)
	{new=(((i<<1)-1)*x*old-(i-1)*older)/(i);
	older=old;old=new;/*forward recursion-stable for all x it seems */
	}
return new;
}


double ql0(l,x) double x;
{
int i;
double log(),old,older,new,qlm();
if(abs(x)>1.) return qlm(l,0,x);/* this sum unstable large x*/
/*sum=0.;twol=l<<1;
for(i=l-1;i>=0;i-=2)
	{
	j=l-i;
	term= (twol-(j))/((double)(l*j));
	sum+= pl0(i,x)*term;
	}
return pl0(l,x)*.5*log(abs((1.+x)/(1.-x)))-sum;*/
older= .5*log(abs((1.+x)/(1.-x)));
if(l==0)return older;
old=older*x-1.;
if(l==1)return old;
for(i=2;i<=l;i++)
	{new=(((i<<1)-1)*x*old-(i-1)*older)/(i);
	older=old;old=new;/*forward recursion-stable for abs(x)<1 ? */
	}
return new;
}

double qlm(l,m,x) double x; int l,m;
{          int i;
double sqrt(),ql0(),pow(),gamma(),f21(),term,log();
if(m==0 && x<1. )return ql0(l,x);
if( abs(x)<1. && l>=0 && m>=0)
	{
	return ((m-1-l)*x*qlm(l,m-1,x)-(m+l-1)*qlm(l-1,m-1,x))/sqrt(1.-x*x);
	/* will ultimately get to ql0 call*/
	}
if(abs(x)==1.)return errorcode;
if(!m)
	{if(l>=0 && l<3)
		{term=.5*log((x+1.)/(x-1.));
		if(!l) return term;
		else if(l==1) return term*x-1.;
		else return term*.5*x*(3.*x-1.)-1.5*x;
		}
	}
/*if(abs(x)>=2.)*/
	{
	term= 1.;
	if(l)for (i=1;i<=l;i++)term*= .5*((i<<1)+1);
	return  pow(2.,(double)(-l))*(m%2?-1.:1.)*
	pow(x*x-1., m*.5)*pow(x,-m-l-1.)*gamma(l+m+1.)/term*
	f21( 1.+.5*(l+m),.5*(1.+l+m),1.5+l, 1/(x*x));
	}
/*else
	{
	term= sqrt(x*x-1.);
	f=f21(.5+m,.5-m,1.5+l, (term-x)/term*.5);
	return f*sqrt(.5*pi)*(m%2?-1.:1.)*
	pow((x-term),.5+l)*gamma(1.+l+m)/(gamma(1.5+l)*sqrt(term));
	}*/
/*return  pow(2.,-1.-l)*sqrt(pi)*(m%2?-1.:1.)*
pow(x*x-1., m*.5)*pow(x,-m-l-1.)*gamma(l+m+1.)/gamma(1.5+l)*
f21( 1.+.5*(l+m),.5*(1.+l+m),1.5+l, 1/(x*x));*/
}



double pli(l,z) int l; double z;
{/*omitting factor of -1^(n/2)= i^n*/
double pow(),sum,zz,term,gamma(),coef,twol;
int m,i,j,k;m=l>>1;twol= l<<1;
if(!l)return 1.;
if(z==0.)
	{if(l%2)return 0.;
	term= ((l/2))%2? -1.:1.;
	for(i=1;i<l;i+=2) term*= i/((double)(i+1));
	return term;

	}
sum=gamma(l+1.);
coef= gamma((l<<1)+1.)/(sum*sum); sum=0.;
zz=z*z;term = 1.;k=1;sum=1.;
for(k=1;k<=m;k++)
	{
	j=k<<1;
	term*= (l-k+1.)*(l-j+2.)*(l-j+1.)/(zz*(k)*(twol-j+2)*(twol-j+1) );
	sum+=term;
	}
return pow(z*.5,(double)l)*sum*coef;
}

#define itmax 50
#define tol 1.e-6

double qli(l,z) int l; double z;
{/*omitting factor of (-i)^n+1       */
double pow(),atan(),sum,zz,term,coef,gamma();int twol,twolp;
int i,j,k;twol= l<<1;twolp=twol+1;
if(!l){if(z==0.)return -pi*.5;return -atan(1/z);}
if(abs(z)<1.e-7)
	{
	if(l%2)
		{if(l==1)return -1.;
		term= ((l+1)/2)%2? -1.:1.;
		for(i=2;i<l;i+=2) term*= i/((double)(i+1));
		return term;
		}
	else{return -pi*.5*pli(l,0.);}

	}
if( abs(z)>1.)
	{
	sum= gamma(l+1.);coef= gamma(twol+2.)/(sum*sum);
	sum=0.;
	zz=z*z;term = 1.;
	for(k=1;k<=itmax;k++)
		{
		j=k<<1;
		sum+=term;
		if( abs(term)<tol*abs(sum))break;
		term*= -(l+k)*(l+j-1)*(l+j)/(zz*(k)*(twolp+j)*(twolp+j-1) );
		}
	return sum/(z*pow(z*.5,(double)l)*coef);
	}
else if(z==1.) return 1.e60;
/*else*/
sum=0.;/* Hochstadt p.158 only terms for which l+k odd*/
j= (l%2)?0:1;term=1.;
for(k=j;k<l;k+=2)
	{
	sum+=term*(((k<<1)+1)<<1)*pli(k,z)/((l-k)*(l+k+1));
	term=-term;
	}
/*printf(" sum=%e\n",sum);*/
return pli(l,z)*(pi*.5-atan(z))-sum;
}

/*static int nsmall;*/

double legendrea(m,n,x,real) int m,n,real;double x;
{/* P m,n real=1 real 0 imaginary arg.*/
int i,j,k; double p,z,w,y,gamma[41],sqrt();
k=n-m;
if(k<0)return 0.;
if(!n)return 1.;
w=z=1.;
if(n!=m)for(i=1;i<=k;i++) z*=x;
gamma[0]=1.;
for(i=1;i<=(n+n);i++)
	{
	gamma[i]=w*gamma[i-1];w+=1.;
	}
if(x==0.)
	{
	i=k>>1;
	if((i<<1)!=k)return 0.;
	p= gamma[m+n]/(gamma[i]*gamma[m+i]);
	}
else
	{
        w=1.;y=1./(x*x+1.e-20);
	if(real){y=-y;w=-w;}
	j=3;p=0.;
	for(i=1;i<=12;i++)
		{
		if( (k+2)/2<i ) break;
		p+=gamma[(n<<1)-(i<<1)+2]*z
			/(gamma[i-1]*gamma[n-i+1]*gamma[k-(i<<1)+j-1]);
		z*=y;
		}
	}
z=1.;
for (i=1;i<=n;i++)z+=z;
p/=z;
if(!real)
	{
	i= n-((n>>2)<<2);/* n-4(n/4) corrected by george from herndon */
	if(i>1)p=-p;
	}
if(!m)return p;
j=m>>1;z=abs(w+x*x);
if(m!=(j<<1))
	{z=sqrt(z);
	j=m;
	}
for(i=1;i<=j;i++)p*=z;
return p;
}

/* for Qn,m (x) x>=1 returns array for n 0 to nmax R ratios of q*/

qleg(m,nmax,x,real,r,q) int m,nmax,real; double x,r[],q[];
{double t,q0,s,log(),sqrt(),atan();int i,n;
if(x==1. && real){fprintf(stderr," x==1 out of range for real qleg\n");
	for(i=0;i<=nmax;i++)q[i]=errorcode;return ;}
if(nmax<=13)n=25;
else n=nmax+7;r[0]=0.;
if(real)
	{
	if(!m) q[0]=.5*log((1.+x)/abs(x-1.));
	else
		{
		t=-1./sqrt(abs(x*x-1.));
		q0=0.;
		q[0]=t;
		for(i=2;i<=m;i++)
			{
			s=(x+x)*(i-1)*t*q[0]+(3.*i-i*i-2)*q0;
			q0=q[0];q[0]=s;
			}
		}
	if(x==1.)q[0]=1.e60;/*big*/
	r[n+1]=x-sqrt(abs(x*x-1.));
	for (i=n;i>=1;i--)
		{
		r[i]=(i+m)/(((i<<1)+1.)*x+(m-i-1)*r[i+1]);
		}
	goto end;
	}
/* imaginary x*/
if(!m)
	{
	if(x<.5) q[0]=atan(x)-.5*pi;
	else q[0]=-atan(1./x);
	}
else
	{
	t=1./sqrt(x*x+1.);
	q0=0.;
	q[0]=t;
	for(i=2;i<=m;i++)
		{
		s= (x+x)*(i-1)*t*q[0]+(3*i-i*i-2)*q0;
		q0=q[0];q[0]=s;
		}
	}
r[n+1]= x-sqrt(x*x+1.);
for(i=n;i>=1;i--)
	{r[i]= (i+m)/((i-m+1.)*r[i+1]-((i<<1)+1)*x);}
/*	printf(" r[%d]=%e\n",i,r[i]);}*/
for(i=1;i<=nmax;i+=2)
	r[i]=-r[i];
end:
for(i=1;i<=nmax;i++)q[i]=q[i-1]*r[i];
/* for x=0, does not do well on odd order terms*/
if(!m && abs(x)<1.e-8 && nmax>0)
	{/*odd l*/
	if(nmax<1)return;
	q[1]=-1.;if(nmax<=1)return;
	t= 1.;
	for(i=2;i<=nmax;i+=2)
		{
		t*= i/((double)(i+1));
		q[i+1]=t*((i/2)%2?1.:-1.);
		}
	return;
	}
return;
}

double qnu(x,nu,real)int real;double x,nu;
{
double y,f21(),sqrt(),gamma(),z,pow(),sin(),cos(),f;
if(x>1.)
	{
	z=1./(x*x);
	if(!real)z=-z;f= f21(1.+.5*nu,.5*(nu+1.),nu+1.5,z);
	/*printf(" f21=%e\n",f);*/
	return sqrt(pi)*gamma(nu+1.)/(gamma(nu+1.5)*pow( 2.*x,nu+1.))
		*f;
	}
else if(x==1. || x==-1.)return errorcode;
if(!real)return errorcode;/* this limitation removable complex ans*/
z= x*x;y=pi*nu*.5;
return sqrt(pi)*(gamma(1.+.5*nu)*cos(y)/gamma(.5*(nu+1.))
	*f21(.5-.5*nu,1.+.5*nu,1.5,z)-gamma(.5*(1.+nu))*sin(y)*.5/gamma(1.+.5*nu)
	*f21(.5*(1.+nu),-.5*nu,.5,z));

}
/* pmunu and qmunu are not intended for direct use and are not
prototyped in  protom.h */

/* CAVEAT: pmunu NOT FULLY GENERAL- USE COMPLEX VERSION cp */
double pmunu(mu,nu,x) double nu,mu,x;
{int n,m;/* -1< x < 2 */
double gamma(),pow(),f21(),z;
n=nu;m=mu;
if( abs(x)<=1. && (double)n==nu && (double)m==mu)return plm(n,m,x);
z=1.-mu;
if( x==1. || (z-((int)z)==0. && z<=0.)) return errorcode;
if(abs(1.-x)<2.) return f21(-nu,nu+1.,1.-mu,.5*(1.-x))*
	pow(abs((x+1.)/(x-1.)),.5*mu)/gamma(1.-mu);
else
	{
	fprintf(stderr," out of range |1-x|<2 p mu nu\n");
	return errorcode;
	}
}

/* CAVEAT: qmunu NOT FULLY GENERAL- USE COMPLEX VERSION cq */
double qmunu(mu,nu,x) double nu,mu,x;
{int n,m; /* returns Qmn, nu(x) divided by exp(i*pi*mu), which
will be complex if mu not an integer. Use complex versions if needed.*/
double gamma(),pow(),f21(),z,y,d1,d2;
n=nu;m=mu;
if(  (double)n==nu && (double)m==mu)return qlm(n,m,x);
z=.5*(1.-x);
/*if(abs(z)<1.)*/
if(mu==m)return errorcode;
d1=1.+nu-mu;
if( (d1==((int)d1)) && d1<=0.)return errorcode;
d2=1.+mu+nu;
if( (d2==((int)d2)) && d2<=0.)return errorcode;
	{
	y=pow(abs((1.+x)/(1.-x)),.5*mu);return
	.5*( f21(-nu,1.+nu,1.+mu,z)/(y*gamma(d1))
	*gamma(-mu)*gamma(d2)
	+gamma(mu)*f21(-nu,1.+nu,1.-mu,z)*y);
	}
/*else
	{
	fprintf(stderr," out of range |1-x|<2 q mu nu x=%le\n",x);
	return errorcode;
	}
*/
}

double pnu(x,nu) double x,nu;
{
double y;
if(abs(x)<1.) return pmunu(0.,nu,x);
y=1./(x*x);
return (gamma(nu+.5)*pow(2.*x,nu)/gamma(nu+1.)
	*f21(.5-.5*nu,-.5*nu,.5-nu,y) +
	gamma(-.5-nu)/gamma(-nu)*pow(2.*x,-nu-1.)
	*f21(.5*nu+1.,.5*(nu+1.),nu+1.5,y)	)/sqrt(pi);
}

