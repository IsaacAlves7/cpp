/* Bessel function tables as a function of (integral) n

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

jn
in
kn
yn
values returned in array double bessela[100];

*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

double bessela[100];

double jn(x,n) int n;double x;
{
int i,top,even,safety=40,k;
double mult,aux,z,y,j1(),j0(),old,older,sum,
	sqrt(),large=1.e10,ilarge;
/* will fill bessela[] with j(x) for n=0,n*/
mult=1.;
if(n<2)
	{
	if(n>0)
		{
		if(x==0.) return 0.;
		bessela[0]=j0(x);
		sum=bessela[1]=j1(x);
		return sum;
		}
	else if (!n)return bessela[0]=j0(x);
	/* else n<0*/
	if (!((-n)%2))mult=-1.;
	}
if(x==0.)return 0.;
if(n>100){printf(" jn:n=%d exceeds size of bessela[]",n);return errorcode;}
aux=2./x;
n=abs(n);
ilarge=1./large;
if( x > (double)n)
	{/*fwd*/
	bessela[0]=older=j0(x);
	bessela[1]=old=j1(x);
	for(i=1;i<n;i++)
		{
		z= i*aux*old-older;
		older=old;
		old=z;
		bessela[i+1]=z;
		}
	return z*mult;
	}
/* backward*/
even =  n+ (int)(sqrt((double)(n*safety)));
top=    (even>>1)<<1;
sum=0.;
even=0;
z=0.;
y=0.;
old=1.;
for(i=top;i>0;i--)
	{
	older=i*aux*old-z;
	if( (i-1)<=n)bessela[i-1]=older;
	z=old;
	old=older;
	if( abs(z)>large)
		{ /*rescale*/
		old*=ilarge;
		z*=ilarge;
		sum*=ilarge;
		y*=ilarge;
		for(k=i-1;k<=n;k++)bessela[k]*=ilarge;
		}
	if(even)sum+=old;
/*if(even) printf(" sum,old %e %e %d\n",sum,old, i); */
	even  ^= 1;
	if(i==n)y=z;
	}
sum= 2.*sum-old;
for(i=0;i<=n;i++)bessela[i]*=(mult/sum);
return y*mult/sum;
}

double in(x,n) int n;double x;
{
int k,even,top,i,safety=40;
double aux,old,older,large=1.e10,ilarge,y,z;
if(!n){bessela[0]=aux=i0(x);return aux;}
if(n==1){bessela[1]=aux=(i1(x));bessela[0]=i0(x);return aux;}
if(x==0.)return 0.;
n=abs(n);
if(n>100){printf(" kn:n=%d exceeds size of bessela[]\n",n);return errorcode;}
aux=2./x;
ilarge=1./large;
/* backward*/
even =  n+ (int)(sqrt((double)(n*safety)));
top=    (even>>1)<<1;
even=0;
z=0.;
y=0.;
old=1.;
for(i=top;i>0;i--)
	{
	older=i*aux*old+z;
	z=old;
	old=older;
	bessela[i-1]=older;
	if( abs(z)>large)
		{
		old*=ilarge;
		z*=ilarge;
		y*=ilarge;
		for(k=i-1;k<=n;k++)bessela[k]*=ilarge;
		}
	if(i==n)y=z;
	}
bessela[0]=aux=i0(x);
aux/=old;
for(k=1;k<=n;k++)bessela[k]*=aux;
return y*aux;
}

double kn(x,n) double x;int n;
{
int k;
double ans,aux,old,older,k0(),k1();
if(x==0.)return errorcode;
if(n==0){bessela[0]=aux=k0(x);return aux;}
else if(n==1){bessela[1]=aux=k1(x);bessela[0]=k0(x);return aux;}
if(n<0)n=-n;
n=abs(n);
if(n>100){printf(" kn:n=%d exceeds size of bessela[]\n",n);return errorcode;}
aux=2./x;
bessela[0]=older=k0(x);
bessela[1]=old=k1(x);
for(k=1;k<n;k++)
	{
	ans=older+k*aux*old;
	older=old;
	old=ans;
	bessela[k+1]=ans;
	}
return ans;
}

double yn(x,n) double x;int n;
{
int k;
double mult,ans,aux,old,older,y0(),y1();
if(x==0.)return errorcode;
if(n==0){bessela[0]=aux=y0(x);return aux;}
else if(n==1){bessela[0]=y0(x);bessela[1]=aux=y1(x);return aux;}
mult=1.;
if(n<0 &&  (-n)%2 )mult=-1.;
n=abs(n);
if(n>100){printf(" yn:n=%d exceeds size of bessela[]\n",n);return errorcode;}
aux=2./x;
bessela[0]=older=y0(x);
bessela[1]=old=y1(x);
for(k=1;k<n;k++)
	{
	ans=k*aux*old-older;
	older=old;
	old=ans;
	bessela[k+1]=ans;
	}
for(k=0;k<=n;k++)bessela[k]*=mult;
return ans*mult;
}

