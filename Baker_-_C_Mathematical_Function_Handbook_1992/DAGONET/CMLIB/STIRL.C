/* Stirling numbers first and second kind

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
#include "cmlib.h"
#include "protom.h"

double stirl2(int n, int m)
{
if (m<0 || n<0 || m>n ) return errorcode;
if(m==n)return 1;
if( !m) return 0;
if(m==1)return 1;
if(m==2) return pow(2., (double)(n-1)) -1;
return m* stirl2(n-1,m)+ stirl2(n-1,m-1);
}

double stirl1(int n, int m)
{
if (m<0 || n<0 || m>n ) return errorcode;
if(m==n)return 1;
if( !m) return 0.;
if(m==1)return  gamma( (double)(n));
return (n-1)* stirl1(n-1,m)+ stirl1(n-1,m-1);
}

double stirlingf(int n, int m)
{
return stirl1(n,m)*((n-m)%2?-1:1);
}