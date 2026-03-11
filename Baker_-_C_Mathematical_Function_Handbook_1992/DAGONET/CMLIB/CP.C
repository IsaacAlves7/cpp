/*  Legendre function P for complex parameters

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#include "complex.h"

#define eps 1.e-8
/* to= i from */
#define CTI(to,from) {to.x=-from.y;to.y=from.x;}

cp(z,mu,nu,ans) struct complex *z,*mu,*nu,*ans;
{
/* need cf21 such that if c neg. integer, returns F/Gamma(c)*/
double sqrt(),gamma(),x,y;
struct complex a,b,c,d,e,f,one,term1,term2,ccos,csin,p,q;
int j,k;
x=z->x;y=z->y;one.x=1.;one.y=0.;
if(abs(y)<eps && abs(x)<=1.)
	{/* on cut  use Bateman 3.4(10) if mu not an integer*/
/*printf(" on cut\n");*/
	if( abs(abs(x)-1.)<eps)
			{/* |x|=1 */
			/*printf(" |x|=1  %le\n",x);*/
			if(cabs(*mu)>eps)/* mu nonzero*/
				{
				/* x=-1 return zero*/
				if( abs(x+1.)<eps)
								{CMPLX(*ans,0.,0.);return 0;}
				/* if mu is an integer return zero*/
				if(abs(mu->y )<eps && abs(mu->x-(int)(mu->x))<eps)
								{CMPLX(*ans,0.,0.);}
				else{CMPLX((*ans),errorcode,0.);}
				return 0;
				}
			/*	mu==0	*/
			if(x==1.){CMPLX((*ans),1.,0.);}
			if(x==-1.)
				{if(cabs(*mu)<eps)
					{CMPLX((*ans),(((int)nu->x)%2?-1.:1.),0.);}
				else{CMPLX((*ans),errorcode,0.);}
				}
			return 0;
			}
/*printf(" on cut but not at end of cut %le\n",x);*/
	CLET(b,one);CADD(b,b,*nu);
	CMPLX(a,0.,0.);CSUB(a,a,*nu);
	CMPLX(f,1.,0.);CSUB(e,f,(*z));CTREAL(f,e,.5);/*f=.5(1-Z)*/
	CMPLX(c,1,0.);CSUB(c,c,*mu);
	/*printf(" cf call\n");printc(&a);printc(&b);printc(&c);printc(&f);*/
	k=cf21(&a,&b,&c,&f,&e);/*printc(&e);printf(" =2F1 return code %d\n",k);*/
	if(cabs(*mu)>eps)
		{
		CMPLX(a,(1.+x)/(1.-x),0.);CTREAL(b,*mu,.5);cpow(&a,&b,&c);
		CMULT(f,c,e);
		}
	else {CLET(f,e);}/* mu=0 */
	/* if mu is a positive integer, use F/Gamma(c) from f21*/
	if(abs(mu->y)<eps && abs(mu->x - ((int)mu->x))<eps && mu->x>0.)
		{ CSET(ans,f);}
	/* otherwise, divide by GAMMA(1-MU)*/	
	else	{CSUB(a,one,*mu);cgamma(&a,&b,&c);CDIV((*ans),f,b);}
	return 0;
	}
/*printf(" off cut\n");*/
CLET(b,one);CADD(b,b,*nu);
CMPLX(a,0.,0.);CSUB(a,a,*nu);
CMPLX(f,1.,0.);CSUB(e,f,(*z));CTREAL(f,e,.5);/*f=.5(1-Z)*/
CMPLX(c,1,0.);CSUB(c,c,*mu);
/*printf(" cf call\n");
printc(&a);printc(&b);printc(&c);printc(&f);*/
cf21(&a,&b,&c,&f,&e);
/*printc(&e);printf(" =2F1 return code=%d\n",k);*/
CLET(b,one);CLET(c,one); CADD(b,b,*z); CSUB(c,c,*z); CDIV(a,b,c);
CTREAL(b,*mu,.5);cpow(&a,&b,&c);CMULT(f,e,c);
if(abs(mu->y)<eps && abs(mu->x - ((int)mu->x))<eps && mu->x>0.)
	{CSET(ans,f);}	/* mu pos integer do not divide by Gamma */
	else {
	CSUB(a,one,*mu);cgamma(&a,&b,&c);CMULT((*ans),f,b);
	}
return 0;
}

