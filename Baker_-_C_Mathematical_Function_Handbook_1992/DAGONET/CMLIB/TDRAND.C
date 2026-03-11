/* test driver for random number generators and distributions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define urand()  randm()

extern long int seed;
extern int naflag,nflag;
extern long int s1,seed2;
extern int s16,s26,s36;
/* TEST CODE BELOW
*/
main(argc,argv) int argc; char **argv;
{
int i,k,l,N,X;
double urand(),normal(),sd,y,x,xs,ys,rk,meant,m2,mu3,mu4,skew,kurt;
double a,b,sum,sq,mean,var,me,ve,p,q,lambda;
double sum1,sum2,sum3,sum4;
double u;
long int outseed;
seed=1;
naflag=1;
nflag=1;
s1=12345;seed2=67890;
s16=12;s26=23;s36=34;
for(i=0;i<=10000;i++)
	{outseed=seed;
	u=randm();
	}
printf(" u=%f seed = %ld \n",u,outseed);
printf(" results of u16:\n");
for(i=0;i<=1000;i++)
	{u=u16();
	if(i%100 ==0)printf(" u=%f \n",u);
	}
printf(" results of u32:\n");
for(i=0;i<=1000;i++)
	{u=u32();
	if(i%100 ==0)printf(" u=%f \n",u);
	}
/* now test distributions */
mean=0.;sd=1.;
while(1)
{printf(" enter count");scanf("%d",&k);
if(k<=0)break;
printf(" kount=%d\n",k);
sum1=sum2=sum3=sum4=0.;
for (i=0;i<k;i++)
	{
	u=normal(mean,sd,&x,&y);
/*	printf(" return normal %f %f\n",x,y);*/
	sum1=sum1+x+y;
	xs=x*x;ys=y*y;
	sum2=sum2+xs+ys;
	sum3=sum3+xs*x+ys*y;
	sum4=sum4+xs*xs+ys*ys;
	}
printf(" tabulating\n");
rk=.5/k;
meant=sum1*rk;
m2=meant*meant;
var=sum2*rk-m2;
mu3=(sum3-3.*sum1*sum2*rk)*rk+2.*meant*m2;
mu4=sum4*rk-3.*m2*m2+6.*m2*sum2*rk-4.*sum3*sum1*rk*rk;
skew=mu3/(var*sqrt(var));
kurt=(mu4/(var*var)-3.)*.5;
/* mean should equal 0. approximately.
   var    "      "   1.     "         (variance=sd*sd)
  skewness           0.
kurtosis             0.                as defined here, which is
CRC Basic Stat. tables definition.  Others define kurtosis without
the fact of .5 and/or without the -3.
*/

printf(" mean=%f,var=%f,skew=%f,kurt=%f\n",meant,var,skew,kurt);
}/* end while*/
mean=0.;sd=1.;
printf(" now fast normal dist mean=0 var=1\n");
while(1)
{printf(" enter count");scanf("%d",&k);
if(k<=0)break;
printf(" kount=%d\n",k);
sum1=sum2=sum3=sum4=0.;
for (i=0;i<k;i++)
	{
	x=na();y=na();
	sum1=sum1+x+y;
	xs=x*x;ys=y*y;
	sum2=sum2+xs+ys;
	sum3=sum3+xs*x+ys*y;
	sum4=sum4+xs*xs+ys*ys;
	}
printf(" tabulating\n");
rk=.5/k;
meant=sum1*rk;
m2=meant*meant;
var=sum2*rk-m2;
mu3=(sum3-3.*sum1*sum2*rk)*rk+2.*meant*m2;
mu4=sum4*rk-3.*m2*m2+6.*m2*sum2*rk-4.*sum3*sum1*rk*rk;
skew=mu3/(var*sqrt(var));
kurt=(mu4/(var*var)-3.)*.5;
/* mean should equal 0. approximately.
   var    "      "   1.     "         (variance=sd*sd)
  skewness           0.
kurtosis             0.                as defined here, which is
CRC Basic Stat. tables definition.  Others define kurtosis without
the fact of .5 and/or without the -3.
*/
printf(" mean=%f,var=%f,skew=%f,kurt=%f\n",meant,var,skew,kurt);
}/* end while*/
printf(" now exponential dist\n");
mean=0.;sd=1.;
while(1)
{printf(" enter count");scanf("%d",&k);
if(k<=0)break;
printf(" kount=%d\n",k);
sum1=sum2=sum3=sum4=0.;
for (i=0;i<k;i++)
	{
	x=ex();y=ex();
/*	printf(" return normal %f %f\n",x,y);*/
	sum1=sum1+x+y;
	xs=x*x;ys=y*y;
	sum2=sum2+xs+ys;
	sum3=sum3+xs*x+ys*y;
	sum4=sum4+xs*xs+ys*ys;
	}
printf(" tabulating\n");
rk=.5/k;
meant=sum1*rk;
m2=meant*meant;
var=sum2*rk-m2;
mu3=(sum3-3.*sum1*sum2*rk)*rk+2.*meant*m2;
mu4=sum4*rk-3.*m2*m2+6.*m2*sum2*rk-4.*sum3*sum1*rk*rk;
skew=mu3/(var*sqrt(var));
kurt=(mu4/(var*var));
/* mean should equal 0. approximately.
   var    "      "   1.     "         (variance=sd*sd)
  skewness           0.
kurtosis             0.                as defined here, which is
CRC Basic Stat. tables definition.  Others define kurtosis without
the fact of .5 and/or without the -3.
*/

printf(" mean=%f,var=%f,skew=%f,kurt=%f\n",meant,var,skew,kurt);
}/* end while*/
mean=0.;sd=1.;
printf(" now cauchy\n");
while(1)
{printf(" enter count");scanf("%d",&k);
if(k<=0)break;
printf(" kount=%d\n",k);
sum1=sum2=sum3=sum4=0.;
for (i=0;i<k;i++)
	{
	x=ca();y=ca();
/*	printf(" return normal %f %f\n",x,y);*/
	sum1=sum1+x+y;
	xs=x*x;ys=y*y;
	sum2=sum2+xs+ys;
	sum3=sum3+xs*x+ys*y;
	sum4=sum4+xs*xs+ys*ys;
	}
printf(" tabulating\n");
rk=.5/k;
meant=sum1*rk;
m2=meant*meant;
var=sum2*rk-m2;
mu3=(sum3-3.*sum1*sum2*rk)*rk+2.*meant*m2;
mu4=sum4*rk-3.*m2*m2+6.*m2*sum2*rk-4.*sum3*sum1*rk*rk;
skew=mu3/(var*sqrt(var));
kurt=(mu4/(var*var));
/* mean should equal 0. approximately.
   var    "      "   1.     "         (variance=sd*sd)
  skewness           0.
kurtosis             0.                as defined here, which is
CRC Basic Stat. tables definition.  Others define kurtosis without
the fact of .5 and/or without the -3.
*/

printf(" mean=%f,var=%f,skew=%f,kurt=%f\n",meant,var,skew,kurt);
}/* end while*/

/*sum=0.;
for(i=0;i<1000;i++){a=random();sum+=a;}
printf(" average random=%le\n",sum/1000.);
sum=0.;
for(i=0;i<1000;i++){a=u16();sum+=a;}
printf(" average u16=%le\n",sum/1000.);*/
while(1)
	{
	printf(" enter int N int k, int m int x for Hypergeom.cumm\n");
	scanf("%d%d%d%d",&N,&k,&l,&X);
	if(N==0)break;
	printf(" N=%d k=%d l=%d X= %d\n",N,k,l,X);
	printf(" =%le\n", hypergeometric_dist(N,l,k,X));
	}
while(1)
	{
	printf(" enter b, int m int x for Erlang.cumm\n");
	scanf("%le%d%d",&p,&k,&X);
	if(k==0)break;
	printf(" =%le\n", erlang_dist(p,k,X));
	}

/*rm=0.;
while(1)
	{printf(" enter neg binomial int x,p ,iter kt\n");
	scanf("%d%le%d",&N,&p,&k);
	if(k<=0 ||p<=0.)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=neg_binomial(N,p);
		printf(" b=%le %le\n",b,randn);
		rm+=randn;
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;rm/=k;
	q=1.-p;	me= N*q/p;ve=me/(p);
	mean=sum/me;
	var= (sq-sum*sum)/(ve);
	printf(" norm mean,var %e %e sum=%e\n",mean,var,sum);
	}
a=0.;
while(1)
	{printf(" enter integer v chisq,iter kt\n");
	scanf("%d%d",&l,&k);
	if(k<=0 ||l<=0)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=chisq(l);
		printf(" b=%e\n",b);
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;
	mean=sum/l;
	var= (sq-sum*sum)/(2*l);
	printf(" normalized mean,var %e %e lambda=%e\n",mean,var,sum);
	}

while(1)
	{
	printf(" enter int X, lambda for Poisson cumm\n");
	scanf("%d%le",&X,&lambda);
	if(lambda==0.)break;
	printf(" =%le\n", poisson_dist(X,lambda));
	}
while(1)
	{
	printf(" enter int N int X, p for binomial cumm\n");
	scanf("%d%d%le",&N,&X,&lambda);
	if(lambda==0.)break;
	printf(" =%le\n", binomial_dist(X,N,lambda));
	}*/
while(1)
	{
	printf(" enter int X int Y, p for neg binomial cumm\n");
	scanf("%d%d%le",&X,&N,&lambda);
	if(lambda==0.)break;
	printf(" =%le\n", neg_binomial_dist(X,N,lambda));
	}

while(1)
	{printf(" enter lognormal mean,sigma,iter kt\n");
	scanf("%le%le%d",&p,&q,&k);
	if(k<=0 ||p<=0.)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=lognormal(p,q);
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;
	a= exp(q*q*.5);me=p*(a);b=a*a;ve=p*p*b*(b-1.);
	mean=sum/me;
	var= (sq-sum*sum)/(ve);
	printf(" norm mean,var %e %e lambda=%e\n",mean,var,sum);
	}

while(1)
	{printf(" enter extreme a,b,kt\n");
	scanf("%le%le%d",&p,&q,&k);
	if(k<=0 ||p<=0.)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=extreme(p,q);
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;
	me=(p)-.57721*q;ve=q*q*pi*pi/6.;
	mean=sum/me;
	var= (sq-sum*sum)/(ve);
	printf(" norm mean,var %e %e lambda=%e\n",mean,var,sum);
	}

while(1)
	{printf(" weibull b,c,kt\n");
	scanf("%le%le%d",&p,&q,&k);
	if(k<=0 ||p<=0.)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=weibull(p,q);
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;
	me= p*gamma((q+1.)/q);ve=p*p*gamma((q+2.)/q) -me*me;
	mean=sum/me;
	var= (sq-sum*sum)/(ve);
	printf(" norm mean,var %e %e lambda=%e\n",mean,var,sum);
	}


while(1)
	{printf(" enter Pareto c,kt\n");
	scanf("%le%d",&p,&k);
	if(k<=0 ||p<=0.)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=pareto(p);
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;
	me=(p)/(p-1.);ve=p/(p-2.)-me*me;
	mean=sum/me;
	var= (sq-sum*sum)/(ve);
	printf(" norm mean,var %e %e lambda=%e\n",mean,var,sum);
	}

/*while(1)
	{printf(" enter binomial n,p ,iter kt\n");
	scanf("%d%le%d",&N,&p,&k);
	if(k<=0 ||p<=0.)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=binomial(N,p);
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;
	q=1.-p;	me= N*p;ve=N*p*(q);
	mean=sum/me;
	var= (sq-sum*sum)/(ve);
	printf(" norm mean,var %e %e lambda=%e\n",mean,var,sum);
	}
*/
/*while(1)
	{printf(" enter geometric p ,iter kt\n");
	scanf("%le%d",&p,&k);
	if(k<=0 ||p<=0.)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=(double)geometric(p);
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;
	q=1.-p;	me= 1./p;ve=q/(p*p);
	mean=sum/me;
	var= (sq-sum*sum)/(ve);
	printf(" norm mean,var %e %e lambda=%e\n",mean,var,sum);
	}
*/
/*
while(1)
	{printf(" enter hypergom int N,X,n,iter kt\n");
	scanf("%d%d%d%d",&N,&X,&nn,&k);
	if(k<=0 ||N<0)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=hypergeometric(N,X,nn);
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;
	me= (double)(nn*X)/((double)N);
	ve=me*(1.-(double)X/(double)N)*(N-nn)/(N-1.);
	mean=sum/me;
	var= (sq-sum*sum)/(ve);
	printf(" norm mean,var %e %e lambda=%e\n",mean,var,sum);
	}
*/
while(1)
	{printf(" enter lambda for Poisson random,iter kt\n");
	scanf("%le%d",&a,&k);
	if(k<=0 ||a<=0.)break;
	sum=sq=0.;
	for(i=0;i<k;i++)
		{b=poisson(a);
		sum+=b;sq+=b*b;
		}
	sum/=k;sq/=k;
	mean=sum/a;
	var= (sq)/(a+a*a);/* NOTvar abt mean, origin*/
	ve= (sq-sum*sum)/(a);/* NOTvar abt mean, origin*/
	printf(" norm mean,var %e %e lambda=%e\n",mean,var,ve);
	}
return 0;
}