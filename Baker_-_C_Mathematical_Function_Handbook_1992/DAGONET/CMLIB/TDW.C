/*  test driver wigner

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define ab(x) ((x)>0. ? (x):0.)
#define min(a,b) ((a)<(b)? (a):(b))
#define max(a,b) ((a)>(b)? (a):(b))
#define DELTA(a,b) ((a)==(b)?1.:0.)
#define tol  3.e-7

void main()
{double j,j1,j2,m,m1,m2; double sqrt();
double a,b,c,d,e,f,g,h,i;
while(1)
{
printf(" for Wigner enter j j1 j2 m m1 m2 j<0 to break\n");
scanf("%le%le%le%le%le%le",&j,&j1,&j2,&m,&m1,&m2);
if(j<0.)break;
printf(" wigner=%e %e\n",wigner(j,j1,j2,m,m1,m2),wigner3j(j,j1,j2,-m,m1,m2));
}
while(1)
{
printf(" for Clebsh j1 m1 j2 m2 j m j<0 to break\n");
scanf("%le%le%le%le%le%le",&j1,&m1,&j2,&m2,&j,&m);
if(j<0.)break;
printf(" CG=%le %le\n",ClebshGordon(j1,m1,j2,m2,j,m),CG(j1,j2,j,m1,m2,m));
if(j==2.)
	{
	h=0.e20;
	if(j==.5 && m==.5 && j2==j1-.5&& m2==-.5-m1)
		h= m1e(j2-m1-.5)*sqrt((j2-m1+.5)/((j2+1.)*(2.*j2+1.)));


	if(m==2.&&j2==j1-1.&& m2==2.-m1)
	h=m1e(j2+m1+1.)*sqrt(20.*(j2-1.+m1)*(j2-m1)*(j2+m1+1.)*(j2+m1+2.)/
	(2.*j2*(2.*j2+1.)*(2.*j2+2.)*(2.*j2+3.)*(2.*j2+4.)));
	else if(m==1.&&j2==j1-1. && m2==1.-m1)
		h=m1e(j2+m1)*2.*(j2+2.-2.*m1)*
		sqrt(5.*(j2+m1)*(m1+j2+1.)/
		(2.*j2*(2.*j2+1.)*(2.*j2+2.)*(2.*j2+3.)*(2.*j2+4.)));
	else if(m==0.&&j2==j1-1.&& m2==-m1)
	h=m1e(j2+m1)*2.*m1*sqrt(30.*abs(j2-m1+1.)
	*(j2+m1+1.)/(2.*j2*(2.*j2+1.)*(2.*j2+2.)*(2.*j2+3.)*(2.*j2+4.)));
	printf(" should be=%le\n",h);
	}
}
while(1)
{printf(" for racah enter a b c d e f a<-100 to break\n");
 scanf("%le%le%le%le%le%le",&a,&b,&c,&d,&e,&f);
 if(a<-100.)break;
printf(" racah %le 6j %le\n",racah(a,b,c,d,e,f),Wigner6j(a,b,e,d,c,f));
printf(" a+c-f=%le\n",a+c-f);
if(e==.5)
	{if((b==a+.5)&&(d==c-.5))
		{printf("W=%le\n",m1e(a+c-f)*
			sqrt(ab((a-c+f+1.)*(f-a+c)
			/((2.*a+1.)*(2.*a+2.)*(2.*c)*(2.*c+1.)))));
		}
	if((b==a+.5)&&(d==c+.5))
		{printf("W=%le\n",m1e(a+c-f)*
			sqrt(ab((a+c+f+2.)*(a+c-f+1.)
			/((2.*a+1.)*(2.*a+2.)*(2.*c+2.)*(2.*c+1.)))));
		}

	}
if(e==1.)
	{if(b==a && c==d)
		{printf("W=%le\n",m1e(a+c+f-1.)*
			(a*(a+1.)+c*(c+1.)-f*(f+1.))
			/sqrt(ab(4.*a*(a+1.)*(2.*a+1.)*c*(c+1.)*(2.*c+1.))));
		}
	if(b==a+1. && d==c-1.)
		{printf("W=%le\n",m1e(a+c-f)*
			sqrt(ab((f-a+c)*(f-a+c-1.)*(a-c+f+2.)*(a-c+f+1.)
			/(4.*(a+1.)*(2.*a+3.)*(2.*a+1.)*c*(2.*c-1.)*(2.*c+1.)))));
		}
	if(b==a+1. && c==d)
		{printf("W=%le\n",m1e(a+c-f)*
			sqrt(ab((a+c+f+2.)*(a-c+f+1.)*(a+c-f+1.)*(c+f-a)
			/(4.*(a+1.)*(2.*a+3.)*(2.*a+1.)*(c)*(2.*c+1.)*(c+1.)))));
		}
	if(b==a+1. && d==c+1.)
		{printf("W=%le\n",m1e(a+c-f)*
			sqrt(ab((a+c+f+3.)*(a+c+f+2.)*(a+c-f+2.)*(a+c-f+1.)
			/(4.*(a+1.)*(2.*a+3.)*(2.*a+1.)*(2.*c+3.)*(2.*c+1.)*(c+1.)))));
		}
	}
if(e==0.)
	{printf("W=%le\n",m1e(a+c-f)*
		 (a==b?1.:0.)*(c==d? 1.:0.)/sqrt((2.*a+1.)*(2.*c+1.)));
	}
if(e==a+b)
	{printf("sb=%le \n",sqrt(ab(gamma(2.*a+1.)*gamma(2.*b+1.)*gamma(a+b+c+d+2.)
	*gamma(a+b+c-d+1.)*gamma(a+b+d-c+1.)*gamma(c+f-a+1.)*gamma(d+f-b+1.)/(
	gamma(2.*(a+b)+2.)*gamma(c+d-a-b+1.)*gamma(a+c-f+1.)*gamma(a+f-c+1.)
	*gamma(a+c+f+2.)*gamma(b+d-f+1.)*gamma(b+f-d+1.)*gamma(b+f+d+2.)
	))));
	}
}

while(1)
	{printf(" enter abcde for Wigner6j ab0/cde a<0 to break\n");
	scanf("%le%le%le%le%le",&a,&b,&c,&d,&e);
	if(a<0.)break;
	printf("6j=%le ab0/cde=%le\n",Wigner6j(a,b,0.,c,d,e),
	m1e(a+e+c)*DELTA(a,b)*DELTA(c,d)/sqrt((2.*a+1.)*(2.*c+1.)));
	g=a*(a+1.)+b*(b+1.)-c*(c+1.);
	printf(" 6j=%le aa1bbc=%le\n",Wigner6j(a,a,1.,b,b,c),
	-g*m1e(a+b+c)/sqrt(a*(2.*a+1.)*(2.*a+2.)*b*(2.*b+1.)*(2.*b+2.)));
	h=(2.*a-1.)*a*4.*(2.*a+1.)*(2.*a+2.)*(2.*a+3.)
	*(2.*b-1.)*b*(2.*b+1.)*(2.*b+2.)*(2.*b+3.);
	if(h>0.)h=m1e(a+b+c)*2.*(3.*g*(g-1.)-4.*a*(a+1.)*b*(b+1.))/sqrt(h);
	else h=0.;
	printf("6j=%le aa2bbc=%le\n",Wigner6j(a,a,2.,b,b,c),h);
	}
while(1)
	{printf(" enter abcdefghi for Wigner9j  a<0 to break\n");
	scanf("%le%le%le%le%le%le %le %le %le",&a,&b,&c,&d,&e,&f,&g,&h,&i);
	if(a<0.)break;
	printf(" 9j=%le %le\n",
		Wigner9j(a,b,c,d,e,f,g,h,i),X(a,b,c,d,e,f,g,h,i));
	if(i==0.)
		{printf("=%le\n",Wigner6j(a,b,c,e,d,g)/sqrt((2.*c+1.)*(2.*g+1.))
		 *DELTA(g,h)*DELTA(f,c)*m1e(b+d+c+g));
		}
	}

}