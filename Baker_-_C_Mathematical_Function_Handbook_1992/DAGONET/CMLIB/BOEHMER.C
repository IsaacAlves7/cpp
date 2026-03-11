/* Boehmer (generalized Fresnel) integral and relatives
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
Boehmer S/C(x,nu) = integral from x to infinity t^(nu-1) sin/cos(t)dt

boehmer(x,nu,type)	type= 0 C 1 for S (type-1,-2 used for Shi,Chi calc)
Si,Ci,Shi,Chi		exponential integral relatives (see A&S)
Fresnel(x,type)		type=0 C 1 for S  Fresnel Integrals
ba(x,nu,type)		used for asymtotic regime.  do not call directly.

Based on Spanier and Oldham with Shi,Chi from Abramowitz & Stegun.
Corrected sign error in S&O 39:12:14
*/
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#define tol 1.e-8

extern int kt; /* term count. if<0, for asymptotic series else series*/
extern double limit;/* series-asymptotic crossover argument for Boehmer*/

/* type = 0  C ;  =1 for S*/

/* aysmptotic formulae for large x*/
double ba(x,nu,typei) int typei;double x,nu;
	{
	double xn,c,s,a,b,term,y,oldterm;int j,type;
	type=typei;
	y=-1./(x*x);
	if(type<0)/* special procedure for Shi,Chi*/
		{y=-y;
		if(type==-1)type=1;
		else type=0;
		c=cosh(x);s=sinh(x);
		}
	else	{c=cos(x);s=sin(x);}
	kt=0;
	for(j=2,a=1.,term=1.;j<1000;j+=2)
		{
		oldterm=term;
		term *= (j-1-nu)*(j-nu)*y;
		if(abs(term)>abs(oldterm))break;
		a+=term;
		kt--;
		if(abs(term)<tol*abs(a))break;
		}
	for(j=3,b=1.,term=1.;j<1000;j+=2)
		{
		oldterm=term;
		term*= (j-1-nu)*(j-nu)*y;
		if(abs(term)>abs(oldterm))break;
		b+=term;
		kt--;
		if(abs(term)<tol*abs(b))break;
		}
	xn=pow(x,nu-2.);
/*printf(" a, b %le %le c s\n",a,b,c,s);*/
	switch(typei)
		{ case 1:  return xn*(x*c*a+(1.-nu)*s*b);/* S ok*/
		  case 0:  return -xn*(x*s*a-(1.-nu)*c*b);/*sign error S&O*/
		 case -1:  return -xn*(x*c*a+(1.-nu)*s*b);/* S ok*/
		  case -2:  return -xn*(x*s*a+(1.-nu)*c*b);/*sign error S&O*/
		  default:
			  fprintf(stderr," badt type to asympt. boehmer\n");
			  return errorcode;
		}
	}

double boehmer( x, nu,typei)int typei;double x,nu;
	{
	double y,x0,sum,term,factor,p,q,gamma(),si(),ci();int j,flag,type;
	type=typei;
	if(limit<=20.) limit=20.;/* default/floor*/
	if(x<0.)return errorcode;
	/* x<0. would require power of negative x */
	/*if(nu>=1.)return errorcode;*/ /* nu>=1 not allowed, but works?*/
	if(nu==0.)flag=1;
	else flag=0;
	if(nu<0.)
		{/* cannot do gamma(nu) nu=0,-1,-2,...*/
		/* for nu=0 use S= .5*pi-Si, C= -Ci(x)*/
		/*recurrences S(x,nu)=(-C(x,nu+1)-x^nu sin(x))/nu
						  C(x,nu)=(S(x,nu+1)-x^nu cos(x))/nu
		may be used for negative nu to get to positive nu */
/*		if( nu==0.)
				{ return type? .5*pi-si(x):ci(x);}*/
		/*else nu<0.*/
		if(x==0.)return errorcode;
		/* pow(0,nu) for nu<0 infinite*/
		if(type)   return (-boehmer(x,nu+1.,0)-pow(x,nu)*sin(x))/nu;
		/*else C*/ return ( boehmer(x,nu+1.,1)-pow(x,nu)*cos(x))/nu;
		}
	/* trap -inf for Ci,Chi(type==0,-2) at x=0, nu=0(flag)*/
	if(flag && ((!type)||(type==-2)) && x==0.)return errorcode;
	if( abs(x) > limit) return ba(x,nu,typei);
	p=nu*pi*.5;
	if(!flag)x0= gamma(nu)*(type? sin(p):cos(p));
	else x0= type? .5*pi : -log(abs(x))- Egamma ;
	y=-x*x;sum=0.;
	if(type<0)/* special procedure for Shi,Chi*/
		{y=-y;
		if(type==-1)type=1;
		else type=0;
		}
	factor= type? 1.+nu:nu;/* first term 1/(1+nu) S or 1/nu C*/
	sum=flag?(type?1.:0.): 1./factor; term=1.;kt=0;
	/*printf(" sum initialized to %le\n",sum);*/
	for(j=1;j<1000;j++)
		{
		p= type? (j<<1)+1:(j<<1);
		term*= y/((p)*(p-1.));
/*printf(" j=%d p=%le,type=%d sum=%le\n",j,p,type,sum);
printf("term=%le term/(p+nu)=%le p=%le\n",term, term/(p+nu),p);*/
		sum+=term/(p+nu);
		kt++;
		if(abs(term)<tol*abs(sum))break;
		}
	if(x>0.){q= pow(x,nu);if(type)q*=x;}
	else q=0.;
/*printf(" [SC](0,nu)=%le sum=%le power=%le\n",x0,sum,q);*/
/*	if(flag && type) return .5*pi-x*(sum);*/
	return x0- q*sum;
	}

double Fresnel(x,type)double x; int type;
	{
	/* corrected from Spanier and Oldham- they
		would have argument x*x, without pi/2 */
	return .5-boehmer(x*x*pi*.5,.5,type)/sqrt(2.*pi);
	}


/*sine and cosine integrals. si(x) is Si(x) not si= Si-pi/2*/

double Si(x) double x;
	{
	return pi*.5-boehmer(x,0.,1);
	}

double Ci(x) double x;
	{return -boehmer(x,0.,0);
	}
double Shi(x) double x;
	{
	return pi*.5-boehmer(x,0.,-1);
	}

double Chi(x) double x;
	{return -boehmer(x,0.,-2);
	}

