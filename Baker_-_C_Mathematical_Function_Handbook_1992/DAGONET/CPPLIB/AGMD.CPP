/* TEST DRIVER FOR:
Elliptic functions of Complex argument and modulus, using
Arithmetic-Geometric Mean (AGM)

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

*/
#include "Cmatrix.hpp"

#define maxterm 100
#define abstol 1.e-10
#define reltol 1.e-5
#define errorcode -1.e60
#define pi 3.141592653589793238462643383279
#define lim 15.
#define max(a,b) ((a)>(b)?(a):(b))
#define infinite_loop for(;;)

extern unsigned _stklen= 54321U;

extern complex csum;// used for E
extern complex Zsum;
//int incomplete;
complex agm( complex& ,complex& ,complex& );
complex F( complex& ,complex& );
complex K(complex& );
complex E(complex& );
complex Z(complex& ,complex& );
complex E(complex& ,complex& );
extern complex jtheta,amplitude;
void scdn( complex& ,complex& ,complex& ,complex& ,complex& );

complex heuman(complex& phi,complex& m);
void theta(complex& v,complex& m,
complex& t1,complex& t2,complex& t3,complex& t4);

void neville(complex& u,complex& m,
complex& ts,complex& tc,complex& td,complex& tn) ;

complex Weierstrass( complex&, complex&, complex&);
complex InvWeier( complex&, complex&, complex&);
complex ElLog( complex&, complex&, complex&);
complex ElExp( complex&, complex&, complex&);

main()
{  complex z,phi,s,c,d,t1,t2,t3,t4;double x,y;
	complex u;

infinite_loop
	{printf(" K&E,etc.:enter complex k");
	scanf("%le%le",&x,&y);z=complex(x,y);
	if(x==0. && y==0.)break;
	K(z).print(" K=","\n");
	E(z).print(" E=","\n");
	printf(" enter complex phi");scanf("%le%le",&x,&y);phi=complex(x,y);
	F(phi,z).print(" F=","----------------------\n");
	printf(" 2nd, incomplete\n");
	E(phi,z).print(" E(phi,k)=","-----------------\n");
	//phi.print(" Z for phi=","\n");
	Z(phi,z).print(" Z=","------------------\n");
	//z.print(" z for heuman lambda=\n","--------------\n");
	//phi.print(" phi for heuman lambda=\n","--------------\n");
	heuman(phi,z).print(" heuman lambda=","\n");
	}
infinite_loop
	{printf(" Jacobi:enter complex m (not k)");
	scanf("%le%le",&x,&y);z=complex(x,y);
	if(x==0. && y==0.)break;
	printf(" enter complex u");scanf("%le%le",&x,&y);phi=complex(x,y);
	phi.print(" phi=","\n");
	scdn(phi,z,s,c,d);
	s.print(" sn=","\n");
	c.print(" cn=","\n");
	d.print(" dn=","\n");
	jtheta.print(" jtheta=","\n");
	amplitude.print(" amplitude=","\n");
	}
infinite_loop
	{printf(" theta functions:enter complex k");
	scanf("%le%le",&x,&y);z=complex(x,y);
	if(x==0. && y==0.)break;
	printf(" enter complex v ");scanf("%le%le",&x,&y);phi=complex(x,y);
	phi.print(" phi=","\n");
	theta(phi,z,t1,t2,t3,t4);
	t1.print(" t1=","\n");
	t2.print(" t2=","\n");
	t3.print(" t3=","\n");
	t4.print(" t4=","\n");
	u=z*z;//m=k^2
	phi.print(" phi=","\n");
	u.print(" m=","\n");
	neville(phi,u,t1,t2,t3,t4);
	t1.print(" neville t1=","\n");
	t2.print(" neville t2=","\n");
	t3.print(" neville t3=","\n");
	t4.print(" neville t4=","\n");
	}

infinite_loop
	{printf(" Weierstrass: enter complex g2 g3 x");
	scanf("%le%le",&x,&y);
	t2=complex(x,y);
	scanf("%le%le",&x,&y);
	t3=complex(x,y);
	scanf("%le%le",&x,&y);
	z=complex(x,y);
	if(x==0. && y==0.)break;
	u=Weierstrass(z,t2,t3);
	u.print("Weierstrass=","\n");
	z=InvWeier(u,t2,t3);
	z.print(" Inverse Weierstrass=","\n");
	//u=ElLog(z,t2,t3);
	//u.print("EllipticLogarithm=","\n");
	//z=ElExp(u,t2,t3);
	//z.print("EllipticExponential=","\n");
	}
infinite_loop
	{
	printf(" EllipticLog/Exp: enter complex a,b,x");
	scanf("%le%le",&x,&y);
	t2=complex(x,y);
	scanf("%le%le",&x,&y);
	t3=complex(x,y);
	scanf("%le%le",&x,&y);
	z=complex(x,y);
	if(x==0. && y==0.)break;
	u=ElLog(z,t2,t3);
	u.print("EllipticLogarithm=","\n");
	z=ElExp(u,t2,t3);
	z.print("EllipticExponential=","\n");
	}
}

