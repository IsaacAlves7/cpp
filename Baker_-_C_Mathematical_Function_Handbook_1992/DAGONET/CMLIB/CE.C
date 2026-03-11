/*  ce E for complex argument 0 <= m <= 1

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

*/
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#define zerotol 1.e-8
#define eps 1.e-6

ce(u,m,ce) double m; struct complex *u,*ce;
{
double wr,wi,s,t,c,d,m1,s1,c1,d1,ep,em,uphi,vphi,fm,fm1,em1;
struct complex z,cs,cc,cd,cim;int ic;
ep=zerotol; em=-ep;CMPLX(cim,0.,-1.);
wr=u->x;wi=u->y;s=abs(wr);t=abs(wi);
if(s<ep && t < ep)
	{ce->x=0.;ce->y=0.;return 0;}
if(wr>ep && wi>ep)ic=1;/* quadrant*/
if(wr<em && wi>ep)ic=2;
if(wr<em && wi<em)ic=3;
if(wr>ep && wi<em)ic=4;
if(t<=ep && wr>=ep)ic=5;/* 5,6 pure real arg*/
if(t<=ep && wr<=em)ic=6;
if(s<=ep && wi>=ep)ic=7;/*pure imaginary arg*/
if(s<=ep && wi<=em)ic=8;
m1=1.-m;
/*printf(" debug ic=%d\n",ic);*/
if(ic<5)
	{
	CMPLX(z,s,t);
	cjef(&z,m,&cs,&cc,&cd);
/*	printc(&cs);printf("= sn(u+iv)\n");*/
	jef(s,m,&s,&c,&d);
	jef(t,m1,&s1,&c1,&d1);
	uphi= asin(s);
	vphi=asin(s1);
	tef(uphi,m,eps,&fm,&em);
	tef(vphi,m1,eps,&fm1,&em1);
	d= s1/c1;
	CMPLX(*ce,em,(t+d1*d-em1));
	CTREAL(cs,cs,m*s*d);
/*	printc(&cs);printf("= sn(u+iv)scaled\n");*/
	CMULT(cd,cim,cs);
/*	printc(&cd);printf("= cd\n");*/
	CADD(*ce,*ce,cd);
/*	printc(ce);printf("= ce\n");*/
	s=ce->x;t=ce->y;
	switch (ic)
		{
		case 2: s=-s;break;
		case 4: t=-t;break;
		case 3:	s=-s;t=-t;break;
		default: break;
		}
	CMPLX((*ce),s,t);
	/*printc(ce);printf("= ce\n");*/
	return 0;
	}

if(ic<7)
	{
	jef(s,m,&s1,&c,&d);
	uphi=asin(s1);
	tef(uphi,m,eps,&fm,&em);
	if(ic==6) em=-em;
	CMPLX(*ce,em,0.);
	return 0;
	}
jef(t,m1,&s,&c,&d);
uphi=asin(s);
tef(uphi,m1,eps,&fm,&em);
em= (t+d*s/c-em);
if(ic==8) em=-em;
CMPLX(*ce,0.,em);
return 0;
}

