/*
Elliptic functions of Complex argument and modulus, using
Arithmetic-Geometric Mean (AGM)

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

K,E : complete integrals of 1st, 2nd kind
F,E : same, of second kind,incomplete
Z Jacobi Zeta
scdn  Jacobian elliptic functions sn,cn,dn
heuman lambda
jacobi theta functions
neville theta functions

for scdn, the modulus m is given; in all other cases k (m=k^2) is given
the angle phi in radians is first argument of the incomplete integrals and Z

*/
#include "Cmatrix.hpp"

#define maxterm 100
#define abstol 1.e-10
#define reltol 1.e-7
#define zerotol 1.e-14
#define tolphi 1.e-6
#define errorcode -1.e60
#define pi 3.141592653589793238462643383279
#define lim 15.
#define max(a,b) ((a)>(b)?(a):(b))
#define infinite_loop for(;;)
#define rabs(x) ((x)<0.? -(x):(x))

complex tangent ( complex& x);
complex tanh(complex& x);
complex sinh(complex& x);
complex asin(complex&x);
complex acos(complex&x);
complex arctan(complex&x);

complex csum;// used for E
complex Zsum;
//int incomplete;

complex agm( complex& a0,complex& b0,complex& c0)
{complex a,b,c,d;double p2=2.;
a=a0;b=b0;csum = c0*c0;
infinite_loop
	{ c=(a-b)*.5;
	//b.print(" b=","\n");
	d=(a+b)*.5;b= (a*b)^.5 ;a=d;
	if ( c.abs() < abstol ) return a;
	// next two lines for computation of E(k) only
	csum += c*c*p2;
	p2*=2.;
	}
}

complex K(complex& k)
{
// k= sin alpha
complex a,b,c;
if( (k-pi*.5).abs() < zerotol ) return -errorcode;
a=1.;c=k; b= ( 1.- k*k) ^ .5 ;
return pi/(2.*agm( a,b,c));
}

complex E(complex& k)
{
// k= sin alpha
complex a,b,c,d;
a=1.;c=k; b= ( 1.- k*k) ^ .5 ;
d=K(k);
return  d*(1.-.5*csum );
}

int cutphipm=1;

complex F( complex& phi0,complex& k)
{complex a,b,c,d,phi,dphi,pr,dphip,m,add,mult(1.,0.),m2(1.,0.);
double p2=2.,dif,newdif;
 int m2p=0,ins;//nq=0,ios=1,i4=0;
 complex or; double reler=1.e20;//oreler=1.e20;
m=k*k;
if(phi0.real()>=0.)
	{
	phi=phi0;
	m2=1.;
	}
else
	{
	phi=-phi0;m2=-1.;
	}
dif=phi.real();
// reduce phi
if(cutphipm)
	{
	ins= (dif/pi);
	dif -= ins*pi;
	if(dif>  .5*pi) {ins++;dif-=pi;}
	if(dif< -.5*pi) {ins--;dif+=pi;}
	if( rabs(dif)> .5*pi)
		fprintf(stderr," F wrong phi reduction %le\n",dif);
	add=2.*ins*K(k);
	}
if( dif<0.)
	{
	dif=-dif;mult=-1.;
	}
phi=complex( dif, phi.imaginary());
// PBLM  if phi0= pi/2. use complete integral
if(phi.abs() < zerotol ) return complex(0.,0.);
if((phi- .5*pi).abs() < zerotol)
		{
		Zsum= 0.;
		return K(k);
		}
// reduce m
// Zsum will not be computed for the given phi and m
// hence cannot allow m reduction for cases where Z or E needed
/*
if(cutphipm)
	{
	a= 1.+m;
	if( (m).abs() > 1.)
		{
		b= 1./k;
		c=k;//c=sqrt(m);
		a= asin( c * phi0.csin() );
		// Zsum will be modified
		return F(a,b)/c;
		}

	if( m.real() <0.)
		{
		a= 1.-m;//(1+m) if arg -m
		b= -m/a; m/(1+m)
		c= 1./(a^.5);
		d=b^.5;
		//Zsum will be modified
		return c*(K(d)-F( pi*.5-phi0,d) );
		}
	}
*/
a=1.; b= ( 1.- m ) ^ .5 ;Zsum= 0.;//c0 irrelevant;
//m.print(" m for F=","\n");
//k.print(" k for F=","\n");phi.print(" phi for F=","\n");
infinite_loop
	{ c=(a-b)*.5;
//b.print(" b=","\n");
	dphi=    arctan(b/a* tangent(phi));
//dphi.print(" dphi=","\n");
//phi.print(" phi=","\n");
dphip=dphi;

//dphip.print("before while dphi=","\n");
while ( dphip.real()<(phi.real()-2.*pi)){m2p++;dphip += 2.*pi;}
//dphip.print("after while dphi=","\n");
if( dphip.real()<(phi.real()-pi)){dphip += pi;}
//dphip.print("after if1 dphi=","\n");
dif= rabs(phi.real()-dphip.real() );
newdif= rabs(phi.real()-dphip.real()-pi );
if( newdif<dif) dphip +=pi;
//dphip.print(" actual dphi used=","\n");
	phi +=dphip;
	or=pr;//oreler=reler;

	d=(a+b)*.5;b= (a*b)^.5 ;a=d;
	reler=(phi-2.*dphip).abs();
	//if ( (c.abs() < abstol) && reler < tolphi ) break;
	if ( (c.abs() < abstol) ) break;
	//for E
	Zsum += c * phi.csin();
//c.print(" c=","\n");
//phi.print(" phi=","\n");
//Zsum.print(" zsum now=","\n");
	p2*=2.;

	}
// Z= Zsum  = E- (E/K)*F .  E= csum+ (E/K) F () are complete
reler=(phi-2.*dphip).abs();
if( reler > 1.e-4)
	fprintf(stderr,"F:poor convergence of phi %le\n",reler);
Zsum= (m2*mult)*Zsum;
return  m2*(phi/(p2*a)*mult+add);
}

complex Z(complex& phi,complex& k)
{complex dum;int savpm;
savpm=cutphipm;
cutphipm=0;
dum= F(phi,k);
cutphipm=savpm;
return Zsum;
}

/* we overload the name E*/

complex E(complex& phi,complex& k)
{complex f,term,total;int savpm;
savpm=cutphipm;
cutphipm=0;
f= F(phi,k);
//f.print(" F again"," ");phi.print(" for phi="," ");
//k.print(" k=","\n");
term=f*E(k)/K(k);
total=Zsum+term;
cutphipm=savpm;
return total;
}

complex jtheta,amplitude;

complex argm( complex& u, complex& m)
{
double bigk,q,p;complex bk; int n;
// would it be safe to do this kind of reduction for other m?
//also, desirable to do split sn(x+iy), etc.
// to avoid u.cexp() blowing up
if( rabs(m.imaginary()) < zerotol && rabs(m.real())<1.)
	{
	bk=(1.-m)^.5;
	bk=K(bk);
	bigk=bk.real()*4.;
	q=u.imaginary();
	n= q/bigk;
	q-=n*bigk;
	p=u.real();
	bk=K(m^.5);
	bigk=bk.real()*4.;
	n= p/bigk;
	p-= n*bigk;
	bk=(1.-m)^.5;
	bk= K( bk);
	if( (u-bk*complex(0.,1.)).abs()< zerotol)
		return errorcode;
	return complex(p,q);
	}
else return u;
}


#define NAGM 40

void scdns( complex& u,complex& m,complex& sn,complex& cn,complex& dn)
{
int i,l,n;
complex a[ NAGM ],c[ NAGM ],phi[ NAGM ],v,am,am1,b,argu,t;
double tol=1.e-7;
/* Jacobian theta function valid for all real m<=1 */
double twon; complex term,k,sum,lou,sqrtm1;
//printf(" entered scnd\n");
//if(m<0.)
//	{/* A&S 16.10.1-4*/
//	t=-m;mu= t/(1.+t);mu1=1./(1.+t);k=sqrt(mu1);v=u/k;
//	jef(v,mu,&term,&sum,&lou);
//	*dn= 1./lou;
//	*sn= k* term* *dn;
//	*cn= sum* *dn;
//	return;
//	}
//if(m>1.)
//	{ /* A&S 16.11.1-4*/
//	mu=1./m; k=sqrt(mu);v=u/k;
//	jef(v,mu,&term,&sum,&lou);
//	*sn=k*term;
//	*cn=lou;
//	*dn=sum;
//	return;
//	}
//v=argm(u,m);
v=u;
if ( v.real()==errorcode)
	{sn=errorcode; cn=sn; dn=sn;
	jtheta=errorcode;
	return ;
	}
am=m;
if(m.abs()==0.)
	{
	sn=v.csin();
	cn=v.ccos();
	dn=complex(1.,0.);
	jtheta=0.;//q=1.
	fprintf(stderr," m==0 return\n");
	return;
	}
else if(m.real()==1. && m.imaginary() == 0.)
	{
	sn=tanh(v);
	cn=1./  sinh(v);
	dn=cn;
	jtheta=0.;
//q=1 using 1+2 sum from 1 to inf -1^k coskz,z=2theta,Jolley 429
	fprintf(stderr," m==1 return\n");
	return;
	}
am1=1.-am;
a[0]=1.;
sqrtm1=(am1)^.5;
b=sqrtm1;
/*c[0]=sqrt(m);not used anywhere*/
for(i=1;i< NAGM;i++)
	{
	a[i]=.5*(a[i-1]+b);
	 if( (a[i]).abs()< zerotol) 
		{
		fprintf(stderr," scdn: a[i] vanishing for i=%d\n",i);
		break;
		}
	c[i]=.5*(a[i-1]-b);
	if((c[i]).abs()<tol)break;
	b=(b*a[i-1])^.5;
	}
n=i;
//printf("scdn: n=%d\n",n);
twon=pow(2., (double)(n));

phi[n]=a[n]*v* twon;
sum=0.;
term= .5/twon;
//term.print(" first term","\n");
for(i=n;i>0;i--)
	{
	if( (a[i]).abs()<zerotol)
		{
		//printf(" scdn a[%d]=0\n",i);
		}
//(a[i]).print(" a=","\n");
//(c[i]).print(" c=","\n");
//(phi[i]).print(" phi=","\n");
//((phi[i]).csin()).print(" csin=","\n");
	argu=c[i]*(phi[i]).csin()/a[i];
//argu.print(" argu=","\n");
	t= asin(argu);
	phi[i-1]=.5*(t+phi[i]);
	sum-= term* ((2.*phi[i-1]-phi[i]).ccos()).clog();
	term*=2.;
//sum.print(" sum=","\n");
	}

argu=phi[0]; amplitude=argu;
//argu.print(" amplitude(argu)=","\n");
sn=argu.csin();
cn=argu.ccos();
lou= (phi[1]-argu).ccos();
//lou.print(" lou=","\n");
if( lou.abs() < zerotol)
	dn= errorcode;
else dn= cn/lou;
//tek(0,m,&k,&term);
k= K( m^.5 );
term= argu.ccos();
//term.print(" term now=","\n");
if( term.abs() < zerotol)
	{
	jtheta=errorcode;
	}
else
	{
	jtheta=sum+.5* (2.*sqrtm1/pi*k * lou/ (term)).clog();
	jtheta=jtheta.cexp();
	}
//jtheta.print(" scdn: returning jtheta","\n");
}

void scdn( complex& u,complex& m,complex& sn,complex& cn,complex& dn)
{
complex s1,c1,d1,s,c,d,x,y,dnom,i(0.,1.),amps;
x= u.real();
y=complex(u.imaginary(),0.);
scdns(x,m,s,c,d);amps=amplitude;
dnom= 1.-m;
scdns(y,dnom,s1,c1,d1);
if( rabs(u.real())<zerotol)amps=i*amplitude;
dnom=1./( c1*c1+m*s*s*s1*s1);
sn= dnom*( s*d1+i*c*d*s1*c1);
cn= dnom*( c*c1-i*s*d*s1*d1);
dn= dnom*( d*c1*d1-i*m*s*c*s1);
// cannot return amps for complex argument! jtheta also not reliable
if( rabs(u.real())>zerotol && rabs(u.imaginary())>zerotol)amps=errorcode;
return ;
}

complex heuman(complex& phi,complex& k)
{
complex m,e,ff,ee,kk,kp,m1; m=k*k;m1=1.-m;kp=m1^.5;
kk=K(k);e=E(k);
ff=F(phi, kp);
ee=E(phi, kp);
return 2./pi*(kk*ee-(kk-e)*ff);
}

void theta(complex& v,complex& m,
complex& t1,complex& t2,complex& t3,complex& t4)
{
complex sn,cn,dn,bigk,k,kp,u,srk,srkp;
k=m^.5;
bigk=K(k);
u=2.*bigk*v/pi;
scdns(u,m,sn,cn,dn);
if (jtheta.real()==errorcode)
	{
	t1=complex(errorcode,errorcode);
	t2=t1;t3=t1;t4=t1;
	return;
	}
kp=(1.-m)^.5;srk=(k)^.5;srkp=(kp)^.5;
t1= sn*srk*jtheta;
t2=cn*srk/srkp*jtheta;
t3=dn/srkp*jtheta;
t4=jtheta;
}

void neville(complex& u,complex& m,
complex& ts,complex& tc,complex& td,complex& tn)
{
complex v,bigk,t1,t2,t3,t4,t10,t20,t30,t40,k;
k=m^.5;
bigk=K(k);
v=u*pi/(2.*bigk);
theta(v,m,t1,t2,t3,t4);
if (t1.real()==errorcode)
	{
	ts=complex(errorcode,errorcode);
	tc=ts;td=ts;tn=ts;
	return;
	}
theta(complex(0.,0.),m,t10,t20,t30,t40);
ts=t1*bigk*2./(pi*t20*t30*t40);
tc=t2/t20;
td=t3/t30;
tn=t4/t40;
}
