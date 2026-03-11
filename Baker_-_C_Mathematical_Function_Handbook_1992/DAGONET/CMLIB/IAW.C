/*  integral of Anger-Weber function

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

based upon algorithm (but not code) of 
by R. F Blackburn and D. R. Wilton
Math Note 57 AFWL Kirtland AFB, NM
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"


#define imax 500
#define tol 1.e-10
#define asympl 15.

/* caveat: making topgam too large will cause overflows!*/
#define topgam 160
static double gammi[topgam],gammh[topgam],gammnh[topgam];


/* call with x==0 to initialize*/
static int cginit=0;

gaminit()
	{ double sp,z; int i;
	sp=sqrt(pi);
	/* gammh[z]= gamma(z+.5)*/
	gammh[0]=sp;z=.5;
	for(i=1;i<topgam;i++)
		{
		gammh[i]=gammh[i-1]*z;
		z+=1.;
		}
	/* gammn[z]= gamma(-z)*/
	gammnh[1]= -2.*sp;z=-.5;
	for(i=2;i<topgam;i++)
		{
		z-=1.;
		gammnh[i]=gammnh[i-1]/z;
		}
	/* gammi[z]= gamma(z)*/
	gammi[0]=errorcode;
	gammi[1]=1. ;z=1.;
	for(i=2;i<topgam;i++)
		{
		gammi[i]=gammi[i-1]*z;
		z+=1.;
		}
	}

double gamtab(double x)
	{
	int i,type; double z,sp;
	if(!cginit) gaminit();
	if(x==0.)
		{
		return errorcode;
		}
	else
		{
		i=(int)(2.*x);
		if( abs(2.*x-(double)i) > tol)
			{fprintf(stderr," arg to gamtab %le not half integer\n",x);
			return gamma(x);
			}
		if(i%2)
			{/* half-integral or error*/
			type=1;
			}	
		else type=0;
		}
	switch (type)
		{
		case 0: /* integer*/
			i=(int)x;
			if(i>=topgam)
				{
				fprintf(stderr," out of gamtab range %le\n",x);
				return gamma(x);
				}
			return gammi[i ];
		case 1: /* half-integer*/
			i= (int)(x-.5);
			if(i>=topgam || -i<=-topgam)
				{
				fprintf(stderr," out of gamtab range %le\n",x);
				return gamma(x);
				}
			if(x>0.)return gammh[i];
			return gammnh[-i];
		}
	}


iaw(m,s,ans) int m; struct complex *s,*ans;
{int i,j;
struct complex sum,x,term;
double gamma(),factor,agam;
if (cabs( *s) > asympl && ( s->y!=0. || s->x >=0.) ) goto asympt;
CMPLX(sum,0.,0.);factor=-1.;
CLET(x,*s);
gamtab(0.);
/*infinite_loop
	{
	printf(" enter x for gamtab (0 to end)\n");
	scanf("%le",&factor);
	if(factor==0.)break;
	printf(" gamtab(%le)=%le\n",factor,gamtab(factor));
	}
*/
for (i=0;i<imax;i++)
	{
	agam=.5*(i-m)+1.;j=agam;
	/*printf(" agam=%e %d\n",agam,i);*/
	
	if( !(j<=0 && ((double)j)==agam) )
		{
		agam=gamtab((agam));
		factor= 1./((i+1)*gamtab(1.+.5*(i+m))*agam);
		if(i%2)factor=-factor;
		CTREAL(term,x,factor);
/*printf(" i=%d agam %le factor %le |term|%le |sum| %le\n"
	,i,agam,factor,cabs(term),cabs(sum));*/
		CADD(sum,sum,term);
		if( cabs(term)<cabs(sum)*tol)break;
		}
	CMULT(term,x,*s);
	CLET(x,term);
	}
switch (m%4)
	{
	case 0 : ans->x=sum.x;ans->y=sum.y;return 0;
	case 1 : ans->x=-sum.y;ans->y=sum.x;return 0;
	case 2 : ans->x=-sum.x;ans->y=-sum.y;return 0;
	default: ans->x=sum.y;ans->y=-sum.x;return 0;
	}

asympt:
CMPLX(x,0.,.5*(m+1.)*pi); CTREAL(sum,*s,-2.);CADD(x,x,sum);
cexp(&x,&sum); CLET(x,*s);
CTREAL(x,x,pi);
ctreal(&x,.5,&term);
CDIV( x,sum,term);CTREAL( x,x, -.5);
CMPLX(sum, 0,.5);CADD(sum,sum,x);
CMPLX( x, s->y  , -(s->x)); clog( &x, &term);
term.x -= digamma( .5*(1.+m));
CTREAL(term,term, (m%2? 0.:2. )/pi*.5);
CADD(*ans,sum,term);
return 0;


}
