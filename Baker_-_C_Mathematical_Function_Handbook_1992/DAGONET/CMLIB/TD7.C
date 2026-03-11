/*  Test Driver Chapter 7
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

extern double limit;
extern int kt,iterp;
FILE *fileid,*ploth,*plotd;
main()
	{double x,y,z,nu,table[10]; struct complex arg,ans;
	int i;double erfc,si,ci,q;
	while(1)
		{printf(" enter x, limit for Shi\n");
		scanf("%le%le",&x,&limit);
		if(limit<=0.)break;
		printf(" Shi Chi %le %le\n",Shi(x),Chi(x));
		}
	printf(" enter limit\n");scanf("%le",&limit);
	while(1)
		{printf("enter x,nu<1000.\n");scanf("%le%le",&x,&nu);
		if(x> 1000.)break;
		printf(" x nu %le %le\n",x,nu);
		printf(" Si %le Ci %le Shi %le Chi %le\n", Si(x),Ci(x),Shi(x),Chi(x));
		printf(" %d series S=%le iter \n",kt,boehmer(x,nu,1));
		printf(" %d series C=%le iter \n",kt,boehmer(x,nu,0));
		printf(" %d series FresnelS=%le iter \n",kt,Fresnel(x,1));
		printf(" %d series FresnelC=%le iter \n",kt,Fresnel(x,0));
		if(x>5.)
			{
			printf(" %d asymp  S=%le iter \n",kt,ba(x,nu,1));
			printf(" %d asymp  C=%le iter \n",kt,ba(x,nu,0));
			}
		}
	printf(" ierfc to compare with A&S Ex 5 p. 305\n");
	printf(" %e %e %e \n", ierfc(1.72,1),ierfc(1.72,2),ierfc(1.72,3));
	ierfctable(1.72,3,table);
	printf(" %e %e %e \n", table[1],table[2],table[3]);
	ierfctable(.5,3,table);
	printf(" %e %e %e \n", table[1],table[2],table[3]);
	fileid=fopen("PLOT.DAT","w");
	fprintf(fileid," 4 \n");
	ploth=fopen("PLOTH.DAT","w");
	fprintf(ploth," 2 \n");
	plotd=fopen("PLOTD.DAT","w");
	fprintf(plotd," 1 \n");
	q=1/.56418958;
	printf(" FWDierfc %e %e\n",ierfcf(.1,2),ierfcf(.5,2));
	printf(" FWDierfc %e %e\n",ierfcf(.1,4),ierfcf(.5,4));
	printf(" FWDierfc %e %e\n",ierfcf(.1,6),ierfcf(.5,6));
	/* proof backward is BAD for small z*/
	printf(" scaled ierfc %e %e\n",q*ierfc(.1,1),q*ierfc(.5,1));
	printf(" scaled ierfcf %e %e\n",q*ierfcf(.1,1),q*ierfcf(.5,1));
	printf(" scaled ierfc %e %e\n",4.*ierfc(.1,2),4.*ierfc(.5,2));
	printf(" scaled ierfcf %e %e\n",4.*ierfcf(.1,2),4.*ierfcf(.5,2));
	printf(" scaled ierfc %e %e\n",32.*ierfc(.1,4),32.*ierfc(.5,4));
	printf(" scaled ierfcf %e %e\n",32.*ierfcf(.1,4),32.*ierfcf(.5,4));
	/* proof forward is BAD for large z*/
	printf(" scaled ierfc %e %e\n",384.*ierfc(.1,6),384.*ierfc(.5,6));
	printf(" scaled ierfcf %e %e\n",384.*ierfcf(.1,6),384.*ierfcf(.5,6));

	printf(" scaled ierfc %e %e\n",122880.*ierfc(.1,10),122880.*ierfc(.5,10));
	printf(" scaled ierfcf %e %e\n",122880.*ierfcf(.1,10),122880.*ierfcf(.5,10));

	printf(" scaled ierfc %e %e\n",32.*ierfc(1.,4),32.*ierfc(5.,4));
	printf(" scaled ierfcf %e %e\n",32.*ierfcf(1.,4),32.*ierfcf(5.,4));

	printf(" scaled ierfc %e %e\n",384.*ierfc(1.,6),384.*ierfc(5.,6));
	printf(" scaled ierfcf %e %e\n",384.*ierfcf(1.,6),384.*ierfcf(5.,6));

	printf(" scaled ierfc %e %e\n",122880.*ierfc(1.,10),122880.*ierfc(5.,10));
	printf(" scaled ierfcf %e %e\n",122880.*ierfcf(1.,10),122880.*ierfcf(5.,10));

	infinite_loop
		{
		printf(" enter x<1000,n for ierfc\n");scanf("%le%d",&x,&i);
		if(x>1000.)break;
		printf(" ierf=%le fwd=%le\n",ierfc(x,i),ierfcf(x,i));
		}

	infinite_loop
		{
		printf(" enter x<1000,y for cerror\n");scanf("%le%le",&x,&y);
		if(x>1000.)break; CMPLX(arg,x,y); cerror(&arg,&ans,1.e-6);
		printc(&ans);printf(" = complex error\n");
		if( y<0.)
			{/* note that w(x+iy) is not conj. of w(x-iy)*/
			arg.y=-arg.y;
			cerror(&arg,&ans,1.e-6);
			q= exp( y*y-x*x)*2.;
			ci= cos(2.*x*y)*q;si=-sin(2.*x*y)*q;
			/* remember to conjugate w(x+iy) and that y >0 here*/
			printf(" should be real=%le imag=%le\n",ci-ans.x,si+ans.y);
			}
		}
	infinite_loop
		{
		printf(" enter x<1000,y for cerfc\n");scanf("%le%le",&x,&y);
		if(x>1000.)break; CMPLX(arg,x,y); cerfc(&arg,&ans);
		printc(&ans);printf(" = complex complementary error\n");
		if(y==0.) printf(" erf=%le\n",1.-ans.x);
		}

	for (i=0;i<50;i++)
		{
		x=(i+1)*.1;
		y=erf(x,&erfc);
		fresnel(x,&ci,&si);
		z=dawson(x);
		printf(" x=%f %f %f %f %f %e %d\n",x,y,erfc,ci,si,z,iterp);
		fprintf(fileid,"%f %le %le %le %le\n",x,y,erfc,ci,si,Chi(x),Shi(x));
		fprintf(ploth,"%f %le %le\n",x,Chi(x),Shi(x));
		fprintf(plotd,"%f %le\n",x,z);
		};
	return (0);
	}
