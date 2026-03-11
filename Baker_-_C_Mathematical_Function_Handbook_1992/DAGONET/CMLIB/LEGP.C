/*  
legptable Pn,m for n=m to ntop. used with sjn by spheroidal wave functions.
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include "cmlib.h"
#include "protom.h"

double leg[100];

legptable(x,ntop,m)int ntop,m;double x;
{/* table of P(x)nm, n=m to ntop*/
int i,twom;
double sqrt(),y,new,term;
y=(1.-x*x);
if(x<-1.||x>1.)y=sqrt(-y);
else y=sqrt(y);
twom=m<<1;
if(!m)new=1.;/* P 0,0*/
else new=-y;/*P 1,1*/
/* first find Pm,m*/
term=3.;
for(i=2;i<=m;i++){new*= -y*term;term+=2.;}
/* new = P m,m old= p m-1,m-1*/
leg[0]=new;
leg[1]=(twom+1)*x*new;/*Pm+1,m */
for(i=2;i<=ntop;i++)
	{ leg[i]=((((m+i)<<1)-1.)*x*leg[i-1]-(i-1+twom)*leg[i-2])/(i);}
return 0;
}

