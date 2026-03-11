/*  test driver Spheroidal Wave Functions and friends

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

extern double diag[100],dsph[100],dsphn[100],sphb[100],leg[100]
	,lmn,csquare;

main()
{int i,m,ntop;float xx;double x,cos(),sphjoin(),c,angular(),zz,yy;
double radial(),rad(),sinh(),cosh(),simb();
 int prolate,odd,nmax,ans,p,o,nn,z;float cc;
 while(1)
{printf(" enter c,prolate, odd, order,nmax\n");
scanf("%e%d%d%d%d",&cc,&p,&odd,&z,&nmax);c=cc;
if(c<0.)break; prolate=p;
solve(c,p,odd,z,nmax);
for(ans=0;ans<nmax;ans++)printf(" eigenv=%e\n",diag[ans]);
printf(" enter i for ith eigenvalue\n");
scanf("%d", &o);lmn=diag[o];
nn=odd;
while(nn<z)nn+=2;
nn+=(o<<1);
printf("  nn=%d\n",nn);
/*printf(" enter n\n");scanf("%d",&nn);printf(" n=%d\n",nn);*/
printf(" lambda=%e n=%d\n",lmn,nn);if(nn==100)break;
p=z;csquare=c;if(!prolate)csquare=-csquare;ans=setd(p,nn);
for(p=0;p<ans;p++)printf(" d=%e %e\n",dsph[p],dsphn[p]);
printf(" join factor first kind:%e N=%e rho=%e\n",
	sphjoin(z,nn,1,ans),nmn(z,nn,ans),rho(z,nn,ans));
printf(" Angular:\n");
for(p=0;p<10;p++){x=p/10.;if(prolate)x=cos(p*pi/18.);
	printf(" %f %d %e\n",x,p,angular(z,nn,1,x,0));}
printf(" Radial, 1st kind:\n");
for(p=0;p<10;p++){x=p/5.;if(prolate)x=1.005+.005*p;
	printf(" %f %d %e\n",x,p,radial(z,nn,1,x));}
printf(" Radial, 2nd kind:\n");
for(p=0;p<10;p++){x=p/5.;if(prolate)x=1.005+.005*p;
	printf(" %f %d %e\n",x,p,radial(z,nn,2,x));}
}
/* spherical bessel check*/
for(m=1;m<=10;m++){
x= .5*m;
zz= cosh(x)*(3./(x*x*x)+1./x)-3./(x*x)*sinh(x);
yy= cosh(x)/x;
c= sinh(x)/x-cosh(x)/(x*x);
printf(" x=%e returned:%e  -1:%e -2:%e -3:%e\n",x,simb(x,5),yy,c,zz);
for(i=0;i<=5;i++)
	printf(" %d %e\n",i,sphb[i]);
}
while(1){
printf(" i:enter n, x\n");
scanf("%d%e",&m,&xx);x=xx;
if(xx<0.||m<0)break;
printf(" %e\n",sinb(x,m));
for(i=0;i<=m;i++)
	printf(" %d %e\n",i,sphb[i]);
}
while(1){
printf(" k:enter n, x\n");
scanf("%d%e",&m,&xx);x=xx;
if(xx<0.||m<0)break;
printf(" %e\n",skn(x,m));
for(i=0;i<=m;i++)
	printf(" %d %e\n",i,sphb[i]);
}
while(1){
printf(" j:enter n, x\n");
scanf("%d%e",&m,&xx);x=xx;
if(xx<0.||m<0)break;
printf(" %e\n",sjn(x,m));
for(i=0;i<=m;i++)
	printf(" %d %e\n",i,sphb[i]);
}
while(1){
printf(" y:enter n, x\n");
scanf("%d%e",&m,&xx);x=xx;
if(xx<0.||m<0)break;
printf(" %e\n",syn(x,m));
for(i=0;i<=m;i++)
	printf(" %d %e\n",i,sphb[i]);
}
while(1){
printf(" P;enter ntop,m x\n");
scanf("%d%d%e",&ntop,&m,&xx);x=xx;
if(xx<0.||m<0)break;
legptable(x,ntop,m);
for(i=0;i<ntop;i++)
	printf(" P(%d,%d)%e\n",m+i-1,m,leg[i]);
}
return 0;
}