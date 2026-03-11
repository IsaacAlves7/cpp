/* real hypergeometric function and relatives.  simple,naive,
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

f21 2F1 (Gauss Hypergeometric)
f12 1F2 
f01 0F1

Note: confluent is 1F1
*/
#include "cmlib.h"
#include "protom.h"

double f21(a,b,c,x)
double a,b,c,x;
{/*hypergeometric function 2F1(a,b;c;x)*/
double sum,term,count,tol=1.e-12;
int top=40,i;
if(x==0.)return 1.;
for(i=0,term=sum=1.;i<top;i++)
	{
	count=(double)(i+1);
	term *=(x*(a+i)*(b+i)/((c+i)*count));
/*	printf(" term=%f sum=%f %d\n",term,sum,i);
*/
	if( abs(term)< tol*abs(sum) )return sum;
	sum+=term;
	}
return(sum);
}

double f12(a,b,c,x)
double a,b,c,x;
{/*hypergeometric function 1F2(a;b,c;x)*/
double sum,term,count,tol=1.e-6;
int top=40,i;
if(x==0.)return 1.;
for(i=0,term=sum=1.;i<top;i++)
	{
	count=(double)(i+1);
	term *=(x*(a+i)/((c+i)*(b+i)*count));
/*	printf(" term=%f sum=%f %d\n",term,sum,i);
*/
	if( abs(term/sum)< tol)return sum;
	sum+=term;
	}
return(sum);
}
double F01(c,x)
double c,x;
{/*hypergeometric function 0F1(;c;x)*/
double sum,term,count,tol=1.e-6;
int top=40,i;
if(x==0.)return 1.;
for(i=0,term=sum=1.;i<top;i++)
	{
	count=(double)(i+1);
	term *=(x/((c+i)*count));
/*	printf(" term=%f sum=%f %d\n",term,sum,i);
*/
	if( abs(term/sum)< tol)return sum;
	sum+=term;
	}
return(sum);
}
