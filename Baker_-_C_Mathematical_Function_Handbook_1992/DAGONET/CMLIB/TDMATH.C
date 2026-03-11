/*  Mathieu functions test driver

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define INDEX(i,j) [(j)+(i)*(coln)]

double cv[120];/* 6 x 20*/

main()
{
int l,n,r,coln,i,j,sol,fnc,norm,k[2],odd;
double q,x,f[3],cvv,fval[3],y;
coln=20;n=5;r=5/*5 odd or 4 even*/;q=5.;
printf(" enter q value\n");scanf("%le",&q);
printf(" q=%le\n",q);
printf(" enter fnc value 1=odd 2=even 3 deriv,odd 4 deriv,even\n");
scanf("%d",&fnc);
printf(" enter r\n");
scanf("%d",&r);
printf(" fnc=%d r=%d\n",fnc,r,n);

/* direct user interface*/
if(fnc%2){/*odd*/n=r;}
else n=r+1;
j=mfcval(n,r,q,cv,coln);
printf(" j=%d\n",j);
for(i=0;i<j;i++)printf("r=n %e\n",cv INDEX(0,i) );

/* a more user friendly interface function: */
odd= fnc%2;
matheign(q,r,odd,cv);
for(i=0;i<j;i++)printf("r=n %e\n",cv INDEX(0,i) );

/* for each order */
for(i=0;i<j;i++)
	{r=i;
	printf(" r=%d\n",i);
	cvv=cv INDEX(0,i);/*q unchanged*/
	sol=3;x= .0*pi;
	printf(" x     norm:neutral  Ince  Stratton\n");
	for(norm=1;norm<4;norm++)
		{
		l=math(x,q,r,cvv,sol,fnc,norm,f,k); fval[norm-1]=f[0];
		}
	printf(" x=%f %le %e %e %d\n",x,fval[0],fval[1],fval[2],l);
	norm=2;
	y=mathieu(x,q,r,cv,sol,fnc,norm);
	printf(" Ince=%le\n",y);

	x=.5*pi;
	for(norm=1;norm<4;norm++)
		{
		l=math(x,q,r,cvv,sol,fnc,norm,f,k); fval[norm-1]=f[0];
		}
	printf(" x=%f %le %e %e %d\n",x,fval[0],fval[1],fval[2],l);
	}
printf(" radial mathieu functions\n");
for(i=0;i<j;i++)
	{r=i;
	cvv=cv INDEX(0,i);/*q unchanged*/
	fnc=2;sol=1;x= .50*pi;/* sol for radial first kind norm irrel*/
	norm=1;
		{
		l=math(x,q,r,cvv,sol,fnc,norm,f,k); fval[norm-1]=f[0];
		}
	printf(" x=%f %le %d\n",x,fval[0],l);
	}
}