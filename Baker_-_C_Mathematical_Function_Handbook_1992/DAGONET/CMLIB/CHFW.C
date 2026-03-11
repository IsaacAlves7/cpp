/* confluent hypergeometric function complex arguments

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

Mwhit	: Whittaker confluent hypergeometric function of the first kind 
Wwhit	:    "            "           "         "      "  "  second kind

*/
#include "cmlib.h"
#include "protom.h"

Mwhit(k,mu,x,ans) struct complex *k,*mu,*x,*ans;
{struct complex a,b,c,d;
CMPLX(a,.5,0.);CADD(a,a,*mu);CSUB(a,a,*k);
CMPLX(b,1.,0.);CTREAL(c,*mu,2.);CADD(b,b,c);
c1f1(&a,&b,x,-1,&d);
if(d.x==errorcode) {ans->x=ans->y=errorcode; return 1;}
CMPLX(c,.5,0.);CADD(c,c,*mu);cpow(x,&c,&a);CMULT(c,d,a);
CTREAL(a,*x,-.5);cexp(&a,&b);CMULT((*ans) ,b,c);return 0;
}

Wwhit(k,mu,x,ans) struct complex *k,*mu,*x,*ans;
{struct complex a,b,c,d;
CMPLX(a,.5,0.);CADD(a,a,*mu);CSUB(a,a,*k);
CMPLX(b,1.,0.);CTREAL(c,*mu,2.);CADD(b,b,c);
cu(&a,&b,x,&d);
if(d.x==errorcode) {ans->x=ans->y=errorcode; return 1;}
CMPLX(c,.5,0.);CADD(c,c,*mu);cpow(x,&c,&a);CMULT(c,d,a);
CTREAL(a,*x,-.5);cexp(&a,&b);CMULT((*ans) ,b,c);return 0;
}
