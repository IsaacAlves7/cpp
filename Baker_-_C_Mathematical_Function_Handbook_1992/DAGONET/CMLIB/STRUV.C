/*  Struve functions general order

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

StruveL L
StruveH H
struve	service routine not intended to be called by user
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define tol 1.e-8
#define tolas 1.e-7
#define xcrit 20.

/* prototype for struve: omit or comment out if not using ANSI compiler*/
double struve( double nu, double x, int l);


double struve(nu,x,l) double nu,x; int l;
	{double y,z,sum,factor,term,oldterm,mult,add;int k;
	if(x<0.)return errorcode;
	if(x>xcrit && abs(nu-((int)nu))<tol)
		{/* attempt to  use asymptotic formula*/
		y=2./x;		
		z=y*y;
		mult=1.;
		if(l){z=-z;mult=-1.; add= in((int)(abs(nu)),z);}
		/* for integer nu, I of order -nu = I of order nu*/
		else add=yn(z, (int)nu);
		for(k=1,sum=term=oldterm=1.;k<1000;k++)
			{
			term *= z*(k-.5)*(.5+nu-k);
			/*printf(" %d oldsum oldterm,term %le %le %le\n",k,sum,oldterm,term);*/
			if( (abs(oldterm)<abs(term) || abs(term)<tol*abs(sum))
				&& k>nu-.5)break;
			sum+=term;
			oldterm=term;
			}
		if( abs(term)>tolas*abs(sum))
			{
			fprintf(stderr,
			" struve: asympt. inaccurate using series term=%le sum=%le\n"
				,term,sum);
			goto series;
			}
		return pow(y,1.-nu)*mult/(sqrt(pi)*gamma(.5+nu))*sum+add;
		}
	series: 
	y=x*.5;
	z=y*y;
	if(!l)z=-z;
	for (k=1,sum=term=1.;k<1000;k++)
		{term *=z/((.5+k)*(.5+k+nu));
		sum+=term;
		if(abs(term)<tol*abs(sum))goto fini;
		}
	fprintf(stderr," struve-series-tolerance not met\n");
	fini:return pow( y, nu+1.)*2./(sqrt(pi)*gamma(nu+1.5))*sum;
	}




double StruveH(nu,x) double nu,x;
	{return struve(nu, x,0);}
double StruveL(nu,x) double nu,x;
	{return struve(nu, x,1);}
