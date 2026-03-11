/*  MacRobert E function

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include "Cmatrix.hpp"

#define maxterm 100
#define abstol 1.e-8
#define reltol 1.e-5
#define errorcode -1.e60
#define pi 3.141592653589793238462643383279
#define lim 15.
#define max(a,b) ((a)>(b)?(a):(b))

complex Gmeijer(int m, int n, int p, int q,Cvector& a,Cvector& c, complex& x);

//MacRobert's E expressed in terms of Miejer's G:
// defining B locally produced a bug 50714 report by Zortech ZTC2B
complex Emacrobert(int p,int q, Cvector& a, Cvector& b, complex& x)
{
Cvector B((1+q));
B.setelement(0,complex(1.,0.));
for(int i=1;i<=q;i++)B.setelement(i, b.element(i-1));
return Gmeijer(p,1,q+1,p,B,a,x);
}
