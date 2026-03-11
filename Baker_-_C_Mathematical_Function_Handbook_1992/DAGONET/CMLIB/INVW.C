/*  Inverse of Weierstrass Elliptic P

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define sig 1.e-6

invp(g2,g3,z,ans)struct complex *z,*ans;double g2,g3;
{double delta,r,q,ss,sqrt(),pp,ppp,lambda,m,acos(),asin(),phi,rt1,rt2,rt3;
int roottp;
struct complex s,s1,s2,r1,r2,r3,factor,num,denom,acoss;
roottp=cubic(0.,-g2*.25,-g3*.25,&rt1,&rt2,&rt3);
if (roottp)
	{
	/*CMPLX(r1,rt1,0.);
	CMPLX(r2,rt2,rt3);
	CMPLX(r3,rt2,-rt3);*/
	/*printf(" delta<0.");*/
	ss= rt1;
	pp= 3.*ss*ss-g2*.25;
	if(pp<0.){fprintf(stderr," pblm p'real root=%e<0\n",ss);return 0;}
	ppp=6.*ss;
	pp=sqrt(pp);
	m= .5-.125*ppp/pp;
	lambda=sqrt(pp);
	/*printf(" lambda=%e m=%e\n",lambda,m);*/
	CLET(num,*z);num.x-=(ss+pp);
	CLET(denom,*z);denom.x+= pp-ss;
	/*printf(" denom=%e %e\n",denom.x,denom.y);*/
	CDIV(s,num,denom);
	/*printf(" calling cacos,x=%e %e\n",s.x,s.y);*/
	cacos(&s,&acoss);
	/*printf(" acoss=%e %e\n",acoss.x,acoss.y);*/
	cef( m,&acoss,&factor,&s,sig);
	CTREAL(*ans,s,.5/lambda);
	return 0;
	}
/* 3 real roots- sort them*/
if(rt1<rt2){q=rt1;rt1=rt2;rt2=q;}
if(rt1<rt3){q=rt1;rt1=rt3;rt3=q;}
if(rt2<rt3){q=rt2;rt2=rt3;rt3=q;}
if( rt1<rt2 || rt1<rt3 ||rt2<rt3){printf(" bad sort %e %e %e\n",rt1,rt2,rt3);
	return 1;}
lambda= rt1-rt3;
if(lambda<0.){fprintf(stderr," invw: bad roots %e %e %e\n",rt1,rt2,rt3);return 1;}
/*printf(" all real roots\n");*/
lambda=.5*sqrt(lambda);
m= (rt2-rt3)/(rt1-rt3);
/*printf(" m= %e lambda=%e\n",m,lambda);*/
CLET(num,*z);num.x-=rt1;
CLET(denom,*z);denom.x-= rt3;
CDIV(s,num,denom);
if(cabs(s)>1.e-10)
	{clog(&s,&factor);CTREAL(factor,factor,.5);cexp(&factor,&s);}
/*printf(" s= %e %e\n",s.x,s.y);*/
cacos(&s,&acoss);
/*printf(" acoss= %e %e\n",acoss.x,acoss.y);*/
cef(m, &acoss,&factor,&s,sig);
/*printf(" cef returned\n");*/
CTREAL(*ans,s,.5/lambda);
return 0;
}



