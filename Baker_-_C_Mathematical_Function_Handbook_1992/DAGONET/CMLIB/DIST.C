/*
Copyright 1991 by Louis Baker. All rights reserved.

a) cumm. prob distributions-
	binomial, negative binomial, Poisson, hypergeometric
b) random numbers for distrib.:
		binomial,neg. binomial,chisq,geometric,hypergeometric,Poisson
		lognormal, weibull,extreme, Pareto

Based largely on Hastings and Peacock, with errors in rng for
 neg. binomial and hypergeometric fixed and chi-sq and
 error hypergeom dist fixed.

*/
#include "cmlib.h"
#include "protom.h"
#include <stdio.h>

#define urand() u16()

/* cummulative binomial, neg. binomial, poisson, hypergeometric dist*/

double binomial_dist(int x,int n, double p)
{
/* also related to incomplete Beta function */
 int i;
 double sum,ratio,q,pow(),term;
 if(p<=0. || p>=1. || x>n || x<0 || n<0 ) return errorcode;
 q=1.-p;
 if(q!=0.)ratio=p/q;
 else ratio=1.;/* q=0 => p=1 */
 term= sum= pow(q,(double)n);
for(i=1;i<=x;i++){term *= (ratio*(n-i+1))/(double)i; sum += term;}
return sum;
}

double neg_binomial_dist(int x,int y, double p)
{
 int i;
 double sum,ratio,q,pow(),term,gamma();
 if(p<=0. || p>=1. || y<0 ||x<1) return errorcode;
 q=1.-p;ratio=q;
 term= sum=1.;
for(i=1;i<=y;i++){term *= (ratio*(x+i-1))/(double)i; sum += term;}
return pow(p,(double)x)*sum;
}

double poisson_dist(int x,double lambda)
{
double sum,exp(),term;
int i;
sum=term=1.;
for(i=1;i<=x;i++){term*= lambda/i ;sum+=term;}
return exp(-lambda)*sum;
}

double hypergeometric_dist(int N, int m, int k, int x)
{double binomial_coef(),sum,term; int r;
double z;
/* N balls in urn, k white; draw m;
cumm prob, i.e. prob that up to x white are drawn H*P has error
in prob dist. function!*/
if(x<0||N<1||k>N||x>m || m>N ) return errorcode;
if(x>k)x=k;/* cummulative prob. includes x>N*/
if(x>m)x=m;
/*printf(" N %d m%d k%d x%d\n",N,m,k,x);*/
/* prob no success f(0)=*/
sum=term= binomial_coef(N-k,m)/binomial_coef(N,m);
if(sum==0.)
	{
	/*fprintf(stderr,
	" hypergeom dist warn:prob of drawing no white balls is zero\n");*/
	for(sum=0.,r=0;r<=x;r++)
		{
		sum+= binomial_coef(k,r)*binomial_coef(N-k,m-r);
		}
	return sum/binomial_coef(N,m);
	}
for(r=0;r<x;r++)
	{
	/*printf(" prob of %d white balls=%le\n",r,term);*/
	z=((double)((r+1)*(N-m-k+r+1)));
	if(z<=0.)
		{/* term does not contribute: m so large, must draw
		more than r white balls. Eg if N=4, k=3, m=2 must draw at least 1 wht
		So	f(0)=0, do not get here.
		*/
		fprintf(stderr,
		" error in hypergeometric distribution N M k x %d %d %d %d\n"
		,N,m,k,x);
		exit(1);
		}
	term*=(m-r)*(k-r)/z;
	sum+=term;
	}
return sum;
}

double erlang_dist(double b,int c,double x)
{double sum,term,factor,exp(); int i;
sum=term=1.;factor=x/b;
for(i=1;i<c;i++)
	{term*= factor/i;
	sum+=term;
	}
return 1.-exp(-factor)*sum;
}


extern double random(),u16(),u32();
extern long int seed,s1,s2; extern int s26,s16,s36;

double randn;

#define smallp  .1
#define nlimbi  0

int binomial(int n, double p)
{
int i,kount;double g,log(),sum,a,q;
q=1.-p;
if(p<=0. || p>=1.) return errorcode;
if(p<smallp && n>nlimbi)
	{
	q= 1./log(1.-p);kount=0;sum=0.;
	while(sum<n)
		{
		a=urand();
		if(a<=0.)
			{/* it is possible that urand=0.but unlikely*/
			continue;}
		g= log(a)*q;sum+=g;kount++;
		}
	}
else
	for(i=kount=0;i<n;i++)if(urand()<p)kount++;
return kount;
}

int neg_binomial(int x, double p)
{/* p success prob.  return number of failures before x successes */
int i,kountl,kountm;double g,log(),sum,a,q;
if(p<=0. || p>=1.) return errorcode;
if(p<smallp && x<nlimbi)
	{
	q= 1./log(1.-p);sum=0.;
	for(i=1;i<=x;i++)
		{
		a=urand(); randn=a;
		if(a<=0.)
			{/* it is possible that urand=0.but unlikely*/
			fprintf(stderr," warn bad urand %le\n",a);
			continue;}
		g= log(a);sum+=g;
		}
	return (int)(sum*q-x +.5);/* round to nearest int*/
	}
else
	{
	kountl=kountm=0;
	while(1)
		{
		a=urand();         
randn=a;
/*printf(" random =%le\n",randn);*/
		if(a<=p)kountl++;/* a<p= failure p=1 all failures*/
		else   kountm++;/*a>=p a sucess */
/*		if(kountm>=x)return kountl;*//* x sucesses-return failure kt*/
		if(kountl>=x)return kountm;  /* x failures-return success????*/
		}
	}
return  ierrorcode;
}

extern nflag;

double chisq(int v)
{
int limit,i;
double add,norm(),urand(),prod=1.,log();
/* call nflag only once for maximum efficiency if will be
many calls to chisq or to norm().  That is, pull stmt from below! */

/* method adding normal^2  works well but is expensive*/
/*
add=0.;
for(i=1;i<=v;i++){prod=norm(0.,1.);randn=prod;add+=prod*prod;}
return add;
*/

/*nflag=1;*/
if( v%2)
	{
	limit=(v-1)>>1;/* -1 not really needed, as shift will do it*/
	add= norm(0.,1.);add*=add;
	}
else
	{
	limit= v>>1;
	add=0.;
	}
for(prod=1.,i=1;i<=limit;i++)
	prod*= urand();
if(prod<=0.) return errorcode;
return add-2.*log(prod);/* H&P give .5 not 2.*/

}

double extreme (double a, double b)
{
double log(),urand(),y,z;
y=0;while(y==0.)y=urand();/*nonzero y*/
z=-log(y);/* as 0<y<=1, log(y) <0 . z>0 */
return a+b*log(z);
}

int geometric(double p)
{double log();return  (int)(log(urand())/log(1.-p)+.99999999);}

int hypergeometric(int N, int X, int smalln)
{/* N balls, X white draw smalln without reeplacement
 return count of white balls drawn */
int i,n,kount=0;double p,d,a;
p=((double)X)/((double)N);/* prob. of drawing a white ball*/
n=N;/* number of balls to draw from*/
for (i=1;i<=smalln;i++)
	{
	a=urand();
	d=(a>=p)?1:0;/* what did we draw? 1=black, 0=white*/
	if(!d)kount++;/* equiv, if a<p*/
/*	p= (n*p-d)*(N-i); WRONG?????*/
	p=( p*n-1.+d)/((double)(n-1));/*new prob. of white ball draw*/
	n--;/*=(1-d);not unchanged. is this function of d,ie n-=d or n-= 1-d?*/
	}
return kount;
}

double lognormal(double median, double sigma)/* sigma=shape parameter*/
{
double exp(),norm();
/* see comments on nflag in chisq() */
/*nflag=1;*/
return median*exp(sigma*norm(0.,1.));
}

double pareto( double c)
{double urand(),pow();
return pow( urand(), -1./c);
}


int poisson(double lambda)
{
int x;
/*  */
double f,urand(),exp(),a,sum;
sum=f= exp(-lambda);a=urand();
for(x=0; x<1000;x++)
	{
	if( a<sum)return x;
	f*= lambda/(x+1);
	sum += sum+f;
	}
return x;
}

double weibull(double b, double c)
{
double urand(),pow(),log(),d;
d=0.;
while(d<=0.)d= urand();
return b*pow(-log(d), 1./c);
}

double binomial_coef(int n, int m)
{ double M,N,a,b,c;
if(m==n)return 1.;
if(!m)return 1.;
if(m>n)return 0.;/*term does not contribute. */
M=m+1;N=n+1;
/*printf(" M,N in b_coef: %le %le\n",M,N);*/
a=gamma(N);b=gamma(M);c=gamma( abs(N-M)+1.);
/*printf(" a,b,c %le %le %le\n",a,b,c);*/
return a/(b*c);
/*return gamma(N)/(gamma(M)*gamma(abs(N-M)));*/
}

