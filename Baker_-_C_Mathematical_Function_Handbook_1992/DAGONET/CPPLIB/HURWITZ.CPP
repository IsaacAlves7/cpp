/*  Riemann & Hurwitz (generalized Riemann) Zeta Function
	series 
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
#include "complex.hpp"
#define pi 3.14159265358979323
#define errorcode -1.e60
#define max(a,b) ((a)>(b)?(a):(b))
#define ABS(x) ((x)<0.?-(x):(x))


prep( complex& u, complex& v, complex& f)
{
f=0.;
do	{
	f+=u^(-v);u++;
	if(v.real()<1.)
		{ while(u.real()>2.){u-=1.;f-=u^(-v);}}
	}while(u.real()<=v.real());
return 0;/* keep compiler happy*/
}

int j,k,kused,jused;

#define tolr 1.e-4
#define tola 1.e-10

complex complex::hurwitz(complex ain)
	{complex v,factor,sum1,sum2,a,b,c,ev,bi,term;
	// sum1,sum2 automatically zeroed to start
	int i;double z,bernoulli(int),scale;
	v=  *this;ev=-v;a=ain;
	prep( a,v,sum1);j=k=4;scale=.25;
	jused=j=max(j, (int)((v.abs()+.999)*scale));
	if(x<0.)jused=j=0;
	b=a+ ((double)j); c=a;
	bi= complex(1.,0.)/(b*b);
	for(i=0;i<j;i++){sum1+= c^ev;c=c+1.;} 
	//sum1.print(" J sum=","\n");
	k= max((1.-x)*8.,k);
	z=2.;
	factor=(v-1.)*v*bi*.5;
	for(i=1;i<=k;i++)
		{
		term= factor*bernoulli(i<<1);sum2+=term;
		factor = factor*(v+z)*(v+z-1.)*bi/((z+2.)*(z+1.));
		if(factor== complex(0.,0.))break;
		if( term.abs() <sum2.abs() * tolr || term.abs()<tola) break;
		z+=2.;kused=k;
		}
	//sum2.print(" k sum=","\n");
	sum1+=(b^ev)*.5;
	ev=1.-v;
	sum2= -(1.+sum2)/ev*(b^ev);
	return sum1+sum2;
	}

complex Riemann_zeta(complex& s)
{return s.hurwitz(complex(1.,0));}


#define lim 15.

complex cgamma( complex& );

complex xi( double t)
{ complex a(.5,t) ;complex b; b=.5*a;//b=.25+it/2
return -.5*(.25+t*t)*((complex(pi,0.))^(-b))
	 *cgamma(b)*Riemann_zeta(a);
}

main() 
{double x,y;
complex a,b(1.,0.),c(0.,1.);
printf(" size complex=%d\n",sizeof(complex));
a=b+c;
a.print("summand=","\n");
printf(" magnitude=%e\n", a.abs());
b=a/c;
b.print("quotient=","\n");
a=b.cexp();
a.print(" exponential=","\n");
b=a.clog();
b.print(" log of a=","\n");
a=complex(1.,0.);b=complex(4.,0.);
a.print(" a ","");b.print(" b ","\n");
b+=1.;b.print(" b after +=1. =","\n");
while(1)
	{
	printf(" enter x,y for s\n");scanf("%le%le",&x,&y);
	b=complex(x,y);
	printf(" enter x,y for a\n");scanf("%le%le",&x,&y);
	a=complex(x,y);// a=1 for Riemann Zeta function
	c=b.hurwitz(a);c.print(" zeta=","");b.print(" for z=","\n");
	printf("kused=%d jused=%d j,k %d %d\n",kused,jused,j,k);
	}
}
