/*  Complementary Error function for complex arguments.

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

cerfc(z,ans) struct complex *z,*ans;
{/*erfc for complex z*/
struct complex i,x,y,p,q;
x.x= - z->y;x.y= z->x;
cerror(&x, &y,1.e-7);
CMULT(p,*z,*z);CTREAL(p,p,-1.);
cexp(&p,&q);CMULT(*ans,q,y);
return 0;
}
