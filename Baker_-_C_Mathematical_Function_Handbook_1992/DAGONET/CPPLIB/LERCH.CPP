/*  Lerch Phi Transcendent (Generalized Riemann/Hurwitz Zeta)

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
#include "complex.hpp"
#define errorcode -1.e60
#define ABS(x) ((x)>0.?(x):-(x))
#define pi 3.141592653589793238462643383279

#define tolabs 1.e-8
#define tolrel 1.e-6
// were 1.e-5,1.e-3 respectively

int size=80,szused;

complex Lerch(complex z,complex a,complex xx)
{
complex f,t,p,q,*h,hold,tpu,ev;int j,k; double diff,tk;
if( size%2)size++;//size must be even integer
h=new complex[size+1] ;
// t=1.;
t= 1./(a^z);
h[0]=1.;k=0;tpu=xx; ev=complex(1.,0.);
while(k<size){
	k++;
//test t[n]=x^n geometric series 1+x+x^2+...=1/(1-x)
//	t*=xx;
	f=ev+a;
	if(f.abs())t= tpu/(f^z);
	else t=0.;
	tpu*=xx;ev+=1.;
	h[k]=h[k-1]+t;
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
hold=h[0];
delete h;
return  hold;
}

main()
{complex u,v,w,x;double a,b,c,d,e,f;
while(1)
	{
	printf(" enter v,u,x");
	    scanf("%le%le%le%le%le%le",&a,&b,&c,&d,&e,&f);
	if(a==0.&&b==0.&&c==0.&&d==0.)break;
		u=complex(c,d);v=complex(a,b);w=complex(e,f);
		x=Lerch(v,u,w);
		//u= 1./(1.-w);
		x.print(" Lerch=","\n");
		//u.print(" 1/(1-x)=","\n");
	printf(" max k used=%d\n",szused);
	};
}

