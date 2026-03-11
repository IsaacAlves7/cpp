/*  Riemann & Hurwitz (generalized Riemann) Zeta Function
	Pade summation 
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
//#include <stream.hpp>    //for Zortech C++
#include <iostream.h> //for TURBO C++
#include <math.h>
#include "complex.hpp"
#define errorcode -1.e60
#define ABS(x) ((x)>0.?(x):-(x))
#define max(a,b) ((a)>(b)?(a):(b))
#define pi 3.141592653589793238462643383279

double Rzeta2( int k)
// Riemann Zeta function of 2*(k+1) returned
{double z[6]={1.64493406684822643647,1.08232323371113819152,
	1.01734306198444913971, 1.00407735619794433938,
	1.00099457512781808534,1.00024608655330804830};
if(k<6)return z[k]; 
// return 1+ 2^(-2k)+...+ 6^(-2k) in simplified form;
if(k>=30)return 1.;
double fk=-((k+1)<<1); double p2=pow(2.,fk);
return 1.+((p2+pow(3.,fk))*(1.+p2)+pow(5.,fk));
}

void prep( complex& u, complex& v, complex& f)
{
f=0.;
do	{
	f+=u^(-v);u++;
	if(v.real()<1.)
		{ while(u.real()>2.){u-=1.;f-=u^(-v);}}
	}while(u.real()<=v.real());
}

#define tolabs 1.e-10
#define tolrel 1.e-5


int size=80,szused,jused;

complex complex::hurwitz(complex u)
{
complex b,c,f,t,p,q,*h,hold,tpu,ev;int i,j,k; double diff,tk,z,scale;
if( *this==1. || ((u-*this+(*this).abs())==0.&& (*this).imaginary()==0.))return complex(errorcode,0.);
prep(u,*this,f);
//u.print("modified a=","");f.print(" modifed sum=","\n");
ev= -(*this);
j=4;scale=.25;
jused=j=max(j, (int)(((*this).abs()+.999)*scale));
if(x<0.)jused=j=0;
b=u+ ((double)j); c=u;
for(i=0;i<j;i++){f+= c^ev;c+=1.;} 
//f.print(" modifed sum=","\n");
if( size%2)size++;//size must be even integer
h=new complex[size+1] ;
t=complex(-2.,0.);
h[0]=1.;k=0; 
tpu=2.*pi*b;
tpu=1./(tpu*tpu);
while(k<size){
	tk=k<<1;
	t*=(*this+tk)*(1.-*this-tk)*tpu;
	z=Rzeta2(k);
	k++;
	h[k]=h[k-1]+t*z;
	q=0.;j=k;
	do	{
		p=h[j-1];
		if(h[j]==p)
			{
			h[j-1]=-errorcode;
			if(p==-errorcode)h[j-1]=q;
			}
		else h[j-1]= q+1./(h[j]-p);
		q=p;j--;
		}while(j);
	if(!(k%2))
		{
		diff=(hold-h[0]).abs();
		if(diff<tolabs || diff<tolrel * h[0].abs())break;
		hold=h[0];
		}
	};
szused=k;
ev= -(*this);
tpu= 1.+ev;
q=-h[0]*(b^tpu)/tpu;
delete h;
return f+.5*(b^ev)+q;
}

main()
{complex u,v,w;double a,b,c,d;

while(1)
	{
	cout << " enter v,u" ;
	    //scanf("%le%le%le%le",&a,&b,&c,&d);
		cin >> a >> b >> c >> d ;
	if(a==0.&&b==0.&&c==0.&&d==0.)break;
		u=complex(c,d);v=complex(a,b);
		w=v.hurwitz(u);
		 w.print(" Hurwitz=","\n");
	cout  << " max k used=" << szused << ", J used: " << jused << "\n" ;
	};
}

