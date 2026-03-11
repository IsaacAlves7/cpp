/*  

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

extern int lim1,itkt;
FILE *osiev,*oclaus,*of123,*odeb,*omu,*olob,*odi;

main()
{double p,q,x,y,z,d,r;int n,j,i;
osiev=fopen("sievert.out","w");
oclaus=fopen("clausen.out","w");
of123=fopen("f123.out","w");
odeb=fopen("debye.out","w");
omu=fopen("mu.out","w");
olob=fopen("lob.out","w");
odi=fopen("dilog.out","w");
fprintf(of123," 4\n");
	for(i=0;i<100;i++)
		{
		x= .1*i;
		fprintf(of123,"%le %le %le %le\n",x, af(x,1),	
		 af(x,2),af(x,3));
		}
fprintf(osiev," 4\n");
	q= pi/180.;
	for(i=0;i<100;i++)
		{
		x= .1*i;
		fprintf(osiev,"%le %le %le %le\n",x, sievert(x,10.*q),	
		 sievert(x,45.*q),	 sievert(x,90.*q));
		}

fprintf(oclaus," 1\n");
	for(i=0;i<100;i++)
		{
		x= pi*i/100.;
		fprintf(oclaus,"%le %le\n",x, clausen(x));
		}
fprintf(odi," 1\n");
	for(i=0;i<100;i++)
		{
		x= 1.*(i)/100.;
		fprintf(odi,"%le %le\n",x, Dlog(x));
		}
fprintf(odeb," 4\n");
	for(i=0;i<100;i++)
		{
		x= .1*i;
		fprintf(odeb,"%le %le %le %le\n",x, debye(x,1),
		 debye(x,2),debye(x,3));
		}
fprintf(olob," 1\n");
	for(i=0;i<100;i++)
		{
		x= pi*.5*i/100.;
		fprintf(olob,"%le %le\n",x, Lob(x));
		}
fprintf(omu," 1\n");
printf(" enter p, q for mu(x,p,q) plots\n");scanf("%le%le",&p,&q);
	for(i=1;i<100;i++)
		{
		x= i/10.;
		fprintf(omu,"%le %le\n",x, mu(x,p,q));
		}
printf(" for z=0 dilog=%le %le\n",dilog(0.), Dlog(0.));
printf(" for z=1 dilog=%le %le\n",dilog(1.), Dlog(1.));


infinite_loop
	{
	printf(" for mu:enter x>=0,alpha,beta,lim1\n");
	scanf("%le%le%le%d",&x,&p,&q,&lim1);
	printf(" x etc %le %le %le %d\n",x,p,q,lim1);
	if(x<0. || lim1 <1)break;
	printf("mu(%le,%le,%le)=%le\n",x,p,q,mu(x,p,q));
	}
infinite_loop
	{printf(" for Clausen integral enter theta degree\n");scanf("%le",&x);
	if(x<=0.)break;
	y= x*pi/180.;
	printf("clausen=%le %le %le\n",x,y,clausen(y) );
	}
printf(" sievert:\n");
for(i=0;i<10;i++)
	{
	x=i*2.;
	for(j=1;j<9;j++)
		{
		if(j<7)d=j*10.;
		else if(j==7)d=75.;
		else d=90.;
		r=d*pi/180.;
		printf(" sievert x=%f theta=%f %e\n",x,d,sievert(x,r));
		}
	}

infinite_loop
	{printf(" for Debye enter x,n\n");scanf("%le%d",&x,&n);
	if(n<=0)break;
	printf(" %le", debye(x,n) );
	}
infinite_loop
	{printf(" for Abramowitz f enter x,n\n");scanf("%le%d",&x,&n);
	if(n<=0)break;
	printf(" %le", af(x,n) );
	}
infinite_loop
		{printf("enter x for Lobachevsky\n");scanf("%le",&x);
		if(x==0.)break;
		printf(" Lob=%le\n",Lob(x));
		}
for (i=1;i<50;i++)
	{ x= i*.02;y=1.-x;printf(" x=%f 1-x=%f Dilog=%le %le\n"
		,x,y,Dlog(x),dilog(y));
	}
infinite_loop
	{
	printf(" enter z, n for polylog\n");
	scanf("%le%d",&x,&n);
	if(x<=0.||n<=0)break;
	z=x;
	printf(" x=%f n=%d polylog=%le \n"
		,x,n,polylog(z,(double)n));
	printf("   iterations used=%d\n",itkt);
	}
return 0;
}

