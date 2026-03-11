
/* "Abramowitz" functions f1, f2, f3 from Ch.27 

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "complex.hpp"
int itkt;

#define pi 3.141592653589793238462643383279
#define errorcode -1.e60
#define ABS(x)  ((x)>0.?(x):-(x))


complex af2s(complex& z)
{
complex x,ln,sp,sq,sum1,sum2;
sp=sqrt(pi);
x=z;
if(x.abs()>0.)ln=x.clog();
else ln=0.;
sq=x*x;
sum1=(-.000491-.0000235*x);
sum1= x*(sq*(-.3225+x*(-.1477
+x*(.03195+x*(.00328+x*sum1))))-1.);
sum1+=sp*.5*(1.+sq);
sum2=.3333333+sq*(.000198*sq-.01667);
/*sum1= .5*sp-x+.5*sp*sq-.3225*cube-.1477*sq*sq
+.03195*cube*sq+.00328*cube*cube-.000491*cube*cube*x
-.0000235*cube*cube*sq;*/

return .5*(ln*x*sq*sum2+sum1);
}

complex af3s(complex& x)
{
complex ln,sp,sq,sum1,sum2;
sp=sqrt(pi);
if(x.abs()>0.)ln=x.clog();
else ln=0.;
sq=x*x;
sum2=.0833+sq*(-.00278+.000025*sq);
/*sum1= 1.-.5*sp*x+.5*sq-.2954*cube+.1014*sq*sq
+.02954*cube*sq-.00578*cube*cube-.00047*cube*cube*x
-.000064*four*four;*/
sum1=(-.00578+x*(-.00047+.000064*x));
sum1= 1.+x*(-.5*sp+x*(.5+x*(-.2954+x*(.1014+x*(
.02954+x*sum1)))));

return .5*(sum1-ln*sq*sq*sum2);
}

complex af1s(complex& x)
{
int k,itmax=400,in; double dk;
complex sp,sq,ln,sum,term,coef,a,b,aa[2],bb[2],t;
sp=sqrt(pi);
sq=x*x;
if(x.abs()>0.)ln=x.clog();
else ln=0.;
sum=0.;
term=1.;
in=0;
aa[0]=0.;aa[1]=-1.;
bb[0]=-sp;bb[1]=1.5*(1.-.57721556649);
a=0.;b=1.;
for(k=0;k<itmax;k++)
	{
	coef= a*ln+b;
	t=coef*term;
	sum+=t;
/*printf(" k in %d %d a,b %e %e s t %e %e c=%e\n",k,in,a,b,sum,t,coef);*/
	if( t.abs()<1.e-8 || t.abs()<1.e-7*(sum.abs()))break;
	term*=x;
	itkt=k;
	if(k>1)
		{
		dk=k;
		coef=1./((dk+1)*(dk-1.)*dk);
		aa[in]= -2.*aa[in]*coef ;
		bb[in]= (-2.*bb[in]-(3.*dk*dk-1.)*aa[in])*coef;
		}
	a=aa[in];b=bb[in];
	in++;
	if(in>1)in=0;
	}
return .5*sum;
}


complex afa(complex& x,int& n)
{
int i,itmax=300; double k;
complex y,z,sum,term,v,vi,a,old,older,t,oldt;
/*if(x<=1.)return -1.;*/
v= 3.*((.5*x)^.666666);vi=1./v;
sum=0.;
term=1.;a=1.;old=0.;t=1.e10;
for(i=0;i<itmax;i++)
	{
	oldt=t;	
	t=term*a;
	if( t.abs()>oldt.abs()) break;
	sum+=t;
	term*=vi;
	older=old;
	old=a;
	if(!i) a=(3.*n*n+3.*n-1.)/12.;
	else
		{
		k=i+1;
		a=(.5*(n-2*k)*(2*k+3-n)*((2*k+3+2*n)*older
			-(12*k*k+36*k-3*n*n-3*n+25)*old))/(12*(k+2));
		}
	if(t.abs()<1.e-8|| t.abs()<1.e-7*(sum.abs()))break;
	itkt=i;
	}
return sqrt(pi/3.)*((v/3.)^(.5*n))*((-v).cexp())*sum;
}

complex af1(complex& x)
{
if(x.abs()<3.)return af1s(x);
return afa(x,1.);
}
complex af2(complex& x)
{
if(x.abs()<2.25)return af2s(x);
return afa(x,2.);
}
complex af3(complex& x)
{
if(x.abs()<3.)return af3s(x);
return afa(x,3.);
}

complex af(complex& x, int n)
{
complex nn,old,older,oldest;
int m;
if(n==1)return af1(x);
if(n==2)return af2(x);
if(n==3)return af3(x);
oldest=af3(x);older=af2(x);old=af1(x);

for(m=4;m<=n;m++)
	{
	nn=.5*( (m-1)*older+x*oldest);
	oldest=older;older=old;old=nn;
	}
return nn;
}

main()
{int i,n;double d;complex u,v,w,z,x;
for (i=0;i<51;i++)
	{d=(i)*(10./50.);
	x=complex(d,0.);
	u=af(x,1);
	v=af(x,2);
	w=af(x,3);
	z=af(x,4);

	printf(" x=%le\n",d);
	u.print(" 1="," ");
	v.print(" 2=","\n");
	w.print(" 3="," ");
	z.print(" 4=","\n");
	}
for (i=0;i<51;i++)
	{d=(i)*(10./50.);
	x=complex(0.,d);
	u=af(x,1);
	v=af(x,2);
	w=af(x,3);
	z=af(x,4);

	printf(" x=%le\n",d);
	u.print(" 1="," ");
	v.print(" 2=","\n");
	w.print(" 3="," ");
	z.print(" 4=","\n");
	}
}