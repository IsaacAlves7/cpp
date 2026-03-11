/*
brent root finder         for real roots of a real function
			root must be bracketed beforehand a<root<b
			function f must change sign: f(a)f(b)<0

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"

#define min(a,b) (((a)<(b))? (a): (b))
#define max(a,b) (((a)>=(b))? (a): (b))

double brent(a,b,eta,t,f)
double a,b,eta,t, (*f) ();
{
double c,d,e,fa,fb,fc,tol,pwr2(),m,p,q,r,s;
fa=(*f)(a);fb=(*f)(b);
init:
	c=a;fc=fa;d=e=b-a;
ext:
	if (abs(fc)<abs(fb) )
		{
		a=b;
		b=c;
		c=a;
		fa=fb;
		fb=fc;
		fc=fa;
		}
tol=2.*eta*abs(b)+t;
m=.5*(c-b);
if(abs(m)>tol &&fb!=0.)
	{
         if(abs(e)<tol || abs(fa)<abs(fb) )
         		{
			d=e=m;
         		}
         else
         	{
                 s=fb/fa;
                 if(a==c)
			{/*linear interp*/
                 	p=(2.*m*s);
                 	q=1.-s;
                 	}
                 else
			{/*inver. quadr. interp*/
                 	q=fa/fc;
                 	r=fb/fc;
                 	p=s*(2.*m*q*(q-r)-(b-a)*(r-1.))  ;
			q=(q-1.)*(r-1.)*(s-1.);
                 	}
                 if(p>0.) q*=-1.;
                 else     p*=-1.;
                 s=e;
		 e=d;
		 if(2.*p< 3.*m*q-abs(tol*q) && p<abs(.5*s*q))
                 		d=p/q;
                 else
				{d=e=m;}

		}/* else abs()<tol...*/
        a=b;
        fa=fb;
        b+=(  abs(d)>tol ? d :(m>0.?tol:-tol)  );
	fb=(*f)(b);
	if((fb>0.)==(fc>0.)) goto init;
goto ext;
       	}
 return(b);
 }
