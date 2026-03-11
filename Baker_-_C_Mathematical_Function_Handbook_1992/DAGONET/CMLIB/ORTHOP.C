/*  orthogonal polynomials

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.


backp	service routine-not intended for direct user call
Pjacobi		Jacobi polynomials
laguerre	Laguerre
Cgegenbauer	Gegenbauer
Tcheby		Chebyshev, 1st kind
Ucheby		Chebyshev, 2nd kind		
Plegendre	Legendre
Hermite		Hermite
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

static double  apass,bpass;
static int odd;

double backp(n,cc,typep,x)double x; int n,cc,typep;
{
int m;   double a,f,b,c;
switch(cc)
	{
	case 1:f=x;break;
	case 2:f=x*x;break;
	case 3: f=1.-x;break;
	default:printf(" base case=%d\n",cc);return errorcode;
	}
for(a=1.,m=n;m>0;m--)
	{
	switch(cc)
		{
		case 1: /* Gen Laguerre*/
			b=n-m+1.;
			c=m*(apass+m);
			break;
		case 2:/* all else*/
			switch(typep)
				{
				case 1:/* gegenbauer*/
					b=2.*(n-m+1.)*(apass+n+m+odd-1.);
					c= m*(((m+odd)<<1)-1);
					break;
				case 2:/* Chebyshev T*/
					b=2.*(n-m+1.)*(n+m+odd-1);
					c= m*(((m+odd)<<1)-1);
					break;
				case 3:/* Chebyshev U (shifted)*/
					b=2.*(n-m+1.)*(n+m+odd);
					c= m*(((m+odd)<<1)-1);
					break;
				case 4: /* legendre*/
					b=(n-m+1)*(((n+m)<<1)+odd-1);
					c= m*(((m+odd)<<1)-1);
					break;
				case 5: /* Hermite*/
					b=2.*(n-m+1);
					c= m*(((m+odd)<<1)-1);
					break;
				}
			break;
		case 3:/*Jacobi P*/
			b=(n-m+1.)*(apass+bpass+n+m);
			c=2.*m*(apass+m);
			break;
		}
	a=1.-b/c*f*a;
	}
return a;
}

double Pjacobi(alpha,beta,n,x) int n; double alpha,beta,x;
{double gamma(),d;
apass=alpha;bpass=beta;
d= gamma(alpha+n+1.)/(gamma(alpha+1.)*gamma((double)n+1.));
return d*backp(n,3,0,x);
}

double laguerre(alpha,n,x) int n; double alpha,x;
{double gamma(),d;
apass=alpha;
d= gamma(alpha+n+1.)/(gamma(alpha+1.)*gamma((double)n+1.));
return d*backp(n,1,0,x);
}

double Cgegenbauer(alpha,n,x) int n; double alpha,x;
{double gamma(),d;
int m;m=n/2;odd=n%2;
apass=alpha;
d=((m%2)?-1.:1.)*pochhammer(alpha,m+odd)/gamma(m+1.);
if(odd)d*=2.*x;
return d*backp(m,2,1,x);
}

double Tcheby(n,x) int n; double x;
{double gamma(),d;
int m;m=n/2;odd=n%2;
d=((m%2)?-1.:1.);
if(odd)d*=(n)*x;
return d*backp(m,2,2,x);
}

double Ucheby(n,x) int n; double x;
{double gamma(),d;
int m;m=n/2;odd=n%2;
d=((m%2)?-1.:1.);
if(odd)d*=2.*(m+1)*x;
return d*backp(m,2,3,x);
}

double Plegendre(n,x) int n; double x;
{double pow(),gamma(),d;
int m;m=n/2;odd=n%2;
d=((m%2)?-1.:1.)*gamma((double)n+1.)
  /(pow(4.,(double)m)*gamma((double)m+1.)*gamma((double)(m+odd+1)));
if(odd)d*=(m+1.)*x;
return d*backp(m,2,4,x);
}

double Hermite(n,x) int n; double x;
{double pow(),gamma(),d;
int m;m=n/2;odd=n%2;
d=((m%2)?-1.:1.)*gamma(n+1.)
  /(gamma((double)m+1.));
if(odd)d*=2.*x;
return d*backp(m,2,5,x);
}
