/*  polynomial root finders
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

solvq	solve quadratic
ctreal	x^n where x is complex and n is real  
cpow	x^n where x, n are complex
ccubic	solve cubic with complex coefficient
cubic	solve cubic with real coefficients
cquartic solve quartic with complex coefficients
quartic  solve quartic with real coefficients
*/


#include <stdio.h>
#include <stdlib.h>
#include "cmlib.h"
#include "complex.h"
#include "protom.h"
/*static csqrt(x,y) struct complex *x,*y;
{struct complex z;clog(x,&z);CTREAL(z,z,.5);cexp(&z,y);}*/

solvq(b,c,ans1,ans2) struct complex *b,*c,*ans1,*ans2;
{
struct complex disc,ca;
CMULT(disc,*b,*b);CTREAL(ca,*c,4.);CSUB(disc,disc,ca);
csqrt(&disc,&disc);CTREAL(disc,disc,.5);
CTREAL(ca,*b,-.5);
CADD(*ans1,ca,disc);CSUB(*ans2,ca,disc);
return 0;}

static struct complex disc;

/* x^y where y real, x complex*/

ctreal(x, y, ans) struct complex *x,*ans; double y;
{struct complex z;
if( cabs((*x))==0.){CMPLX((*ans),0.,0.);return 1;}
clog(x,&z);CTREAL(z,z,y);cexp(&z,ans);
return 0;
}

cpow(x, y, ans) struct complex *x,*y,*ans;
{struct complex z,a;
if( cabs((*x))==0.){CMPLX((*ans),0.,0.);return 1;}
clog(x,&z);CMULT(a,z,*y);cexp(&a,ans);
return 0;
}

/*  could make eps an external (global) variable for control*/
#define eps 1.e-8
double polyresid;/*residual error for return if of interest*/

ccubic(a1,a2,a3,r1,r2,r3) struct complex *a1,*a2,*a3,*r1,*r2,*r3;
{struct complex h,g,offset,d1,d2,d3,d4,u,v,w,w2,b1,b2,b3,ut,vt;
double third=1./3.,bestr;
int j,k;
CMPLX(w, -.5, sqrt(3.)*.5);
CONJG(w2,w);
CTREAL(offset,*a1,third);
CMULT(d2,*a1,*a1);CTREAL(d1,*a2,3.);CSUB(d1,d1,d2);
CTREAL(h,d1,third*third);
CMULT(d3,d2,*a1);CTREAL(d3,d3,2./27.);
CMULT(d2,*a1,*a2);CTREAL(d2,d2,third);
CSUB(g,d3,d2);
bestr=cabs(h);
CADD(g,g,(*a3));
CMULT(disc,g,g);CTREAL(g,g,-1.);
CMULT(d1,h,h);CMULT(d2,d1,h);CTREAL(d1,d2,4.);/* d2=h^3 d1=4h^3*/
CADD(disc,disc,d1);
/*printc(&disc);printf("=disc\n");*/
csqrt( &disc, &d1);
/*printc(&d1);printf("=sqrt disc\n");*/
CADD(d4,d1,g);CTREAL(u,d4,.5);
/*CSUB(d3,g,d1);CTREAL(v,d3,.5);*/
if(bestr==0.){CMPLX(u,0.,0.);}
CSUB(v,g,u);
/*printc(&h);printf("=h\n");printc(&g);printf("=g\n");
printc(&u);printf("=u\n");printc(&v);printf("=v\n");*/
/*CDIV(v,d2,u);CTREAL(v,v,-1.); pblm if |u| small */
CMULT(d3,u,v);CADD(d3,d3,d2);
/*printc(&d3);printf(" = uv+H^3 should be zero\n");*/
 bestr=1.e37;
 ctreal(&v,third,&d1);CLET(ut,d1);
 ctreal(&u,third,&d1);CLET(vt,d1);
 for(j=0;j<3;j++)
	{
	for(k=0;k<3;k++)
		{
/*		printc(&ut);printf("=u^third\n");
		printc(&vt);printf("=v^third\n");*/
		CADD( *r1,ut,vt);
		CMULT(d1,ut,w);CMULT(d2,vt,w2);
		CADD(*r2,d1,d2);
		CMULT(d1,ut,w2);CMULT(d2,vt,w);
		CADD(*r3,d1,d2);
		CSUB(*r1,*r1,offset);
		CSUB(*r2,*r2,offset);
		CSUB(*r3,*r3,offset);
		/* check*/
		CADD(h,*a1,*r1);CMULT(g,h,*r1);CADD(h,g,*a2);
		CMULT(g,h,*r1);CADD(g,g,*a3);
		polyresid=cabs(g);
/*		printc(r1);printf("=r1, resid %le \n",polyresid);*/
		if( polyresid < eps)return 0;
		if(polyresid<bestr)
			{bestr=polyresid;
			CLET(b1,*r1);CLET(b2,*r2);CLET(b3,*r3);
			}
		CMULT(*r1,vt,w);CLET(vt,*r1);
		}
	CMULT(*r1,ut,w);CLET(ut,*r1);
	}
		{
		/*if(bestr>eps)*/
		fprintf(stderr," ccubic: best resid=%le > eps %le\n"
			,bestr,eps);
		CSET(r1,b1);CSET(r2,b2); CSET(r3,b3);
		polyresid=bestr;
		if(bestr>eps)return 1;
		}
return 0;
}

int cubic(a,b,c,r1,r2,r3) double a,b,c,*r1,*r2,*r3;
{/* solve cubic x^3+ax^2+b^x+c=0.
if 3 real roots, return 0, roots in r1,2,3
if 1 real, two complex. real in r1, r2 is real pt +/-r3*i return1*/
struct complex a1,a2,a3,cr1,cr2,cr3;
CMPLX(a1,a,0.);
CMPLX(a2,b,0.);
CMPLX(a3,c,0.);
ccubic(&a1,&a2,&a3,&cr1,&cr2,&cr3);
if(  abs(disc.y)> eps)
	{fprintf(stderr," complex cubic discriminant: %le %le\n",disc.x,disc.y);
	return ierrorcode;
	}
if( disc.x <=0)
	{ /* if =0, some roots equal but still all are real*/
	*r1= cr1.x;
	*r2= cr2.x;
	*r3= cr3.x;
	return 0;
	}
/* one real, two complex conj roots*/
if( abs(cr1.y ) < eps)
	{ *r1= cr1.x;
	*r2=cr2.x;
	*r3=cr2.y;
	}
else if( abs(cr2.y ) < eps)
	{ *r1= cr2.x;
	*r2=cr1.x;
	*r3=cr1.y;
	}
else 
	{ *r1= cr3.x;
	*r2=cr2.x;
	*r3=cr2.y;
	}
return 1;
}

cquartic(a1,a2,a3,a4,r1,r2,r3,r4)
struct complex *a1,*a2,*a3,*a4,*r1,*r2,*r3,*r4;
{
struct complex a,b,c,d,q1,q2,q3,m,q,g,h,b1,b2,b3,b4;
int try,sel;double bestr;
CTREAL(a,*a2,-.5);
CMULT(b,*a1,*a3);CTREAL(b,b,.25);CSUB(b,b,*a4);
CMULT(c,*a2,*a4);CTREAL(c,c,4.);CMULT(d,*a3,*a3);CSUB(c,c,d);
CMULT(d,*a1,*a1);CMULT(q1,d,*a4);CSUB(c,c,q1);CTREAL(c,c,.125);
try=ccubic(&a,&b,&c,&q1,&q2,&q3);
if(try){fprintf(stderr," cquartic bad return ccubic\n");return 1;}
/*
printf(" q1 %le %le\n",q1.x,q1.y);
printf(" q2 %le %le\n",q2.x,q2.y);
printf(" q3 %le %le\n",q3.x,q3.y);
 */
try=sel=0;bestr=1.e30;
CLET(q,q1);
if( abs(q2.y)/(cabs(q2)+eps) < abs(q.y)/(cabs(q)+eps)){CLET(q,q2);sel=1;}
if( abs(q3.y)/(cabs(q3)+eps) <abs(q.y)/(cabs(q)+eps)){CLET(q,q3);sel=2;}

attempt:
/*printf(" sel=%d try=%d\n",sel,try);*/
CMULT(a,*a1,*a1);CTREAL(a,a,.25);CTREAL(b,q,2.);CADD(a,a,b);
CSUB(a,a,*a2); csqrt(&a,&m);polyresid=cabs(m);
if( polyresid< eps)
	{CMULT(a,q,q);CSUB(a,a,*a4);csqrt(&a,&b);}
else
	{CMULT(b, q,*a1);CSUB(c,b,*a3);CDIV(b,c,m);CTREAL(b,b,.5);}
/* solve quadratics x^2 +(a1*.5 -m)x +q1 -b  =0
		    x^2 +(a1*.5 +m)x +q1 +b =0 */
CTREAL(a,*a1,.5);CSUB(a,a,m);
CSUB(c,q,b);
solvq(&a,&c,r1,r2);
CTREAL(a,*a1,.5);CADD(a,a,m);
CADD(c,q,b);
solvq(&a,&c,r3,r4);
/* check */
CADD(h,*a1,*r1);CMULT(g,h,*r1);CADD(h,g,*a2);CMULT(g,h,*r1)CADD(g,g,*a3);
CMULT(h,g,*r1);CADD(g,h,*a4);polyresid=cabs(g);
if( polyresid > eps)
	{try++;
	if( polyresid <bestr)
		{bestr=polyresid;CLET(b1,*r1);CLET(b2,*r2);CLET(b3,*r3);CLET(b4,*r4);
		}
	/*printf(" retry as resid=%le\n",polyresid);*/
	if(!sel){CLET(q,q2);sel=1;}
	else{CLET(q,q3);sel=2;}
	if(try<3)goto attempt;
	polyresid=bestr;
	CSET(r1,b1); CSET(r2,b2);	CSET(r3,b3); CSET(r4,b4);
	fprintf(stderr," cquartic: resid=%le > teps=%le\n"
		,polyresid,eps);
	return 1;
	}
return 0;
}

quartic(a,b,c,d,r1,r2,r3,r4) double a,b,c,d; struct complex *r1,*r2,*r3,*r4;
{
struct complex a1,a2,a3,a4;int retval;
CMPLX(a1,a,0.);
CMPLX(a2,b,0.);
CMPLX(a3,c,0.);
CMPLX(a4,d,0.);
retval=cquartic(&a1,&a2,&a3,&a4,r1,r2,r3,r4);
/* could use return value to signal type of roots if retval=0*/
return retval;
}
