/*  Struve functions lowest order

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

h0
h1
l0
l1
h0a,l0a,l1a,h1a asymptotic
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

double h0a(x) double x;
{
int i,itmax=30,k;
double sq,sum,rx,term,oldterm;
sum=0.;
rx=1./x;
term=1.;
sq=rx*rx;
for(i=0;i<itmax;i++)
	{
	sum+=term;
	oldterm=term;
	k= (i<<1) +1;
	k *=k;
	term*= -sq*(k);
	if (abs(term)>abs(oldterm) || abs(term)<1.e-6) break;
	}
return 2./pi*rx*sum+y0(x);
}

double h1a(x) double x;
{
int i,itmax=30,k;
double sq,sum,rx,term,oldterm;
sum=0.;
rx=1./x;
term=1.;
sq=rx*rx;
for(i=0;i<itmax;i++)
	{
	sum+=term;
	oldterm=term;
	if(!i) k=-1;
	else
		{
		k= (i<<1)-1;
		k *=(k+2);
		}
	term*= -sq*(k);
	if (abs(term)>abs(oldterm) || abs(term)<1.e-6) break;
	}
return 2./pi*sum+y1(x);
}

double l0a(x) double x;
{
int i,itmax=30,k;
double sq,sum,rx,term,oldterm,i0();
sum=0.;
rx=1./x;
term=-1.;
sq=rx*rx;
for(i=0;i<itmax;i++)
	{
	sum+=term;
	oldterm=term;
	k= (i<<1)+1;
	k *=k;
	term*= -sq*(k);
	if (abs(term)>abs(oldterm) || abs(term)<1.e-6) break;
	}
return 2./pi*rx*sum+i0(x);
}

double l1a(x) double x;
{
int i,itmax=30,k;
double sq,sum,rx,term,oldterm,i1();
sum=0.;
rx=1./x;
term=-1.;
sq=rx*rx;
for(i=0;i<itmax;i++)
	{
	sum+=term;
	oldterm=term;
	if(!i) k=-1;
	else
		{
		k= (i<<1)-1;
		k *=(k+2);
		}
	term*= -sq*(k);
	if (abs(term)>abs(oldterm) || abs(term)<1.e-6) break;
	}
return 2./pi*sum+i1(x);
}

double l0(x)
double x;
{
double term,count,sq,sum,xa=10.,l0a(),tol=1.e-6;
int i,top=20;
if(x==0.)return 0.;
if(x>xa) return l0a(x);
sq=x*x;
for(i=0,sum=0.,count=3.,term=x;i<top;i++)
	{
	sum +=term;
	if( sum)
		{
		if( (i>3) & (abs(term/sum)<tol) )return 2./pi*sum;
		}
	term *=   sq/(count*count);
	count += 2.;
	}
return 2./pi*sum;
}

double l1(x)
double x;
{
double term,count,sq,sum,xa=10.,l1a(),tol=1.e-6;
int i,top=20;
if(x==0.)return 0.;
if(x>xa) return l1a(x);
sq=x*x;
for(i=0,sum=0.,count=3.,term=sq/3.;i<top;i++)
	{
	sum +=term;
	if (sum)
		{
		if((i>3) & (abs(term/sum)<tol) )return 2./pi*sum;
		}
	term *=   sq/(count*(count+2.));
	count += 2.;
	}
return 2./pi*sum;
}
double h0(x)
double x;
{
double term,count,sq,sum,h0a(),xa=10.,tol=1.e-6;
int i,top=30;
if(x==0.)return 0.;
if(x>xa) return h0a(x);
sq=x*x;
for(i=0,sum=0.,count=3.,term=x;i<top;i++)
	{
	sum +=term;
	if( sum)
		{
		if( (i>3) & (abs(term/sum)<tol) )return 2./pi*sum;
		}
	term *=   -sq/(count*count);
	count += 2.;
	}
return 2./pi*sum;
}

double h1(x)
double x;
{
double term,count,sq,sum,xa=10.,h1a(),tol=1.e-6;
int i,top=20;
if(x==0.)return 0.;
if(x>xa) return h1a(x);
sq=x*x;
for(i=0,sum=0.,count=3.,term=sq/3.;i<top;i++)
	{
	sum +=term;
	if (sum)
		{
		if((i>3) & (abs(term/sum)<tol) )return 2./pi*sum;
		}
	term *=   -sq/(count*(count+2.));
	count += 2.;
	}
return 2./pi*sum;
}
