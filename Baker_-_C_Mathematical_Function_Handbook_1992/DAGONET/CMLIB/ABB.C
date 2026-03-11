
/*  Bessel functions and integrals thereof
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
jt	J(x) for order  1/3
jmt	J(x) for order -1/3
it	I(x) for order  1/3
imt	I(x) for order -1/3
jbes	J(x,nu)
jas	J asymptotic large x
ybes	Y(x,nu)
kbes	K(x,nu)
ibes  I(x,nu)
Ai
Bi 
dai
dbi
IAi
IBi
*/
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"
#include "complex.h"

double jt(z) double z;
{
double pow(),sqrt(),x,y,zz,nu,sum,term;
int k,n,itmax=30;
y= pow(z*.5,.33333333);
sum=0.;
zz=z*z*3.;
term=1.;
for(k=0;k<itmax;k++)
	{
	sum+=term;
	term*= -.25*zz/(double)((k+1)*(3*k+4));
	if( abs(term)<1.e-6 || abs(term/sum)<1.e-5)break;
	}
return sum*y/.892979;

}
double jmt(z) double z;
{
double pow(),sqrt(),x,y,zz,nu,sum,term;
int k,n,itmax=30;
y= pow(z*.5,-.33333333);
sum=0.;
zz=z*z*3.;
term=1.;
for(k=0;k<itmax;k++)
	{
	sum+=term;
	term*= -.25*zz/(double)((k+1)*(3*k+2));
	if( abs(term)<1.e-6 || abs(term/sum)<1.e-5)break;
	}
return sum*y/1.354118;
}

double it(z) double z;
{
double pow(),sqrt(),x,y,zz,nu,sum,term;
int k,n,itmax=30;
y= pow(z*.5,.33333333);
x=z*sqrt(3.);
sum=0.;
zz=x*x*.25;
term=1.;
for(k=0;k<itmax;k++)
	{
	sum+=term;
	term*= zz/(double)((k+1)*(3*k+4));
	if( abs(term)<5.e-7 || abs(term/sum)<1.e-6)break;
	}
return sum*y/.892979;

}
double imt(z) double z;
{
double pow(),sqrt(),x,y,zz,nu,sum,term;
int k,n,itmax=30;
y= pow(z*.5,-.33333333);
x=z*sqrt(3.);
sum=0.;
zz=x*x;
term=1.;
for(k=0;k<itmax;k++)
	{
	sum+=term;
	term*= .25*zz/(double)((k+1)*(3*k+2));
	if( abs(term)<5.e-7 || abs(term/sum)<1.e-7)break;
	}
return sum*y/1.354118;
}

static double asympt=6.;

double jbes(z,nu) double z,nu;
{
double pow(),sqrt(),jas(),gamma(),x,y,zz,sum,term,denom;
int k,n,itmax=50;
if( abs(z)>= asympt)return jas(z,nu);
if( abs(nu)<1. && z!=0.) return jbes(z,nu+1.)*2.*(nu+1.)/z-jbes(z,nu+2.);
x=z*.5;
y=1.;
if(nu!=0.)y= pow(x,nu);
sum=0.;
zz=x*x;
term=1.;
for(k=0;k<itmax;k++)
	{denom=gamma(nu+1.+k);
	if( abs(denom)<1.e-6)
		{
		fprintf(stderr," warn denom=%e nu %e k %d in jbes()\n",
			denom,nu,k);
		break;
		}
	sum+=term/denom;
	term*= -zz/( ((double) (k+1)) );
	if( abs(term)<5.e-7 || abs(term)<1.e-7*abs(sum))break;
	}
return sum*y;
}

double ybessel;

double jas(z,nu) double z,nu;
{
double cos(),sin(),sum,sum2,term2,term,x,y,q,arg,c,s;
double pow(),sqrt(),gamma(),zz,f1,f2,f1o,f2o,nuph,numh,npho,mu;
int k,twok,n,itmax=30;
/* INCREASING NU ONLY MAKES IT WORSE*/
/*if( abs(nu)<1.) return jas(z,nu+1.)*2.*(nu+1.)/z-jas(z,nu+2.);*/
y= sqrt(2./(pi*z));
sum=0.;
sum2=0.;
zz=4.*z*z;
term2=.5/z;
term=1.;
nuph=nu+.5;numh=nu-.5;npho=nu+1.5;
f1=1.e10;f2=1.e10;
mu=4.* nu*nu;
/*
for(k=0;k<itmax;k++)
	{
	f1o=f1;
	twok=k<<1;
	q=(double)(twok);
	if(term!=0.)
		{
		f1=term *gamma(nuph+q)/gamma(nuph-q);
		if( abs(f1/f1o)>1.)
			{
			term=0.;
			f1=0.;
			}
		else term*= -1./(zz*((twok+2)*(twok+1)) );
		sum+=f1;
		}
	f2o=f2;
	if(term2!=0.)
		{
		f2=term2 *gamma(npho+q)/gamma(numh-q);
		if(abs(f2/f2o)>1.)
			{term2=0.;
			f2=0.;
			}
		 else term2*= -1./(zz*((twok+3)*(twok+2)) );
		sum2+=f2;
		}
	if( (abs(term)<1.e-7 || abs(term/sum)<1.e-6)&&
	 ( abs(term2)<1.e-7 || abs(term2/sum2)<1.e-6)  )break;
	}
arg= z-pi*(nu*.5+.25);
return y*(cos(arg)*sum-sin(arg)*sum2);
*/
/*HANKEL EXPANSION:*/
zz=1./(64.*z*z);
term=(mu-1.)*(mu-9.)*zz;
sum=1.-term*.5+term*(mu-25.)*(mu-49.)*zz/24.;
term2=(mu-1.)*.125/z;
sum2=term2*(1.-(mu-9.)*(mu-25.)*zz/6. );
arg= z-pi*(nu*.5+.25);
c=cos(arg);s=sin(arg);
ybessel=y*(s*sum+c*sum2);
return y*(c*sum-s*sum2);
}

double ybes(z,nu) double z,nu;
{
int n,iv;
struct complex zz,jj,yy,dummy;
double x;
n=nu;
if( ( nu-((double)(n))) ==0.)
	{/* call bessel return y*/
	zz.x=z;zz.y=0.;
	bessel(n,&zz,&jj,&yy,&dummy,&dummy,&dummy,&dummy,&iv);
	return yy.x;
	}
if( abs(z)>asympt)
	{
	x=jbes(z,nu);
	return ybessel;
	}
x=(jbes(z,nu)*cos(pi*nu)-jbes(z,-nu))/sin(pi*nu);
/*printf(" y returning %e\n",x);*/
return  x;
}

double kbes(z,nuu) double z,nuu;
{
int n;
struct complex zz,ans;
double tol=1.e-6,x,y,sum,term,factor,pow,old,mu,nu;
/* not valid z<0*/
nu=abs(nuu);
if(z<0.)return errorcode;
n=nu;
if( (nu-((double)(n)))==0.)
	{/* call kbess return k*/
	zz.x=z;zz.y=0.;
	kbess(&zz,&ans,n);
	return ans.x;
	}
/* for z<0 fix */
if( (z)>=asympt)
	{
	x= 1./(8.*z);
	mu=nu*nu*4.;
	sum=1.+(mu-1.)*x*((mu-9.)*x*.5*(1.+(mu-25.)*x*.33333)+1.);
	return exp(-z)*sqrt(pi/(2.*z))*sum;
	}
/* below not used in asymptotic regime where sign of nu irrelevant*/
return pi*.5*( ibes(z,-nu)-ibes(z,nu))/sin(nu*pi);
}

/* NOT STATIC- MAKE c1,c2 PUBLIC and usable
c1a = Ai(0)
c2a= -Ai'(0)
*/
double c1a=.355028053887817,c2a=.258819403792807;

double smallf(z) double z;
{
double x,y,zz,nu,sum,term;
int k,n,itmax=50;
sum=0.;
zz=z*z*z;
term=1.;
for(k=0;k<itmax;k++)
	{
	sum+=term;
	y= zz/((double)((3*k+2)*(3*k+3)));
	term*=y;
	if( abs(term)<5.e-7 || abs(term/sum)<5.e-7)break;
	}
return sum;

}
double smallg(z) double z;
{
double pow(),sqrt(),x,y,zz,nu,sum,term;
int k,n,itmax=50;
sum=0.;
zz=z*z*z;
term=1.;
for(k=0;k<itmax;k++)
	{
	sum+=term;
	term*= zz/((double)((3*k+3)*(3*k+4)));
	if( abs(term)<5.e-7 || abs(term/sum)<5.e-7)break;
	}
return z*sum;

}

double ai(z) double(z);
{
int i,itmax=40;
double x,zeta,zi,pow(),sqrt(),arg,exp(),c,d,term,term2,old,old2,
	tol=1.e-6,factor,sum,sin(),cos(),sum2;
if(z==0.) return c1a;
x=abs(z);
if(x>=asympt)
	{
	zeta= x*sqrt(x)*.6666666;
	zi=1./zeta;
	if(z>0.)
		{
		sum=0.;
		c=term=1.;
		old=1.;
		for(i=0;i<itmax;i++)
			{
			sum+=term;
			if(!i)factor=15.;
			else factor*= (6*i-5)*(6*i-3)*(6*i-1)/(2*i-1);
			c*=factor/(216.*(i+1));
			term*=-zi*c;
			if( abs(c)<tol || abs(c)<tol*abs(sum)) break;
			if( abs(term)>abs(old))break;
			old=term;
			}
		return .5/sqrt(pi*sqrt(x))*exp(-zeta)*sum;
/*printf(" asymp Ai %e\n", .5/sqrt(pi*sqrt(x))*exp(-zeta)*sum);*/
		}
else
	{
	arg=zeta+pi*.25;
	sum=sum2=0.;
	d=c=term=1.;term2=-1.;
	old=old2=1.;
	for(i=0;i<itmax;i++)
		{
		if(i%2)sum2+=term2;
		else sum+=term;
		if(!i)factor=15.;
		else factor*= (6*i-5)*(6*i-3)*(6*i-1)/(2*i-1);
		c*=factor/(216.*(i+1));
		term*=-zi*c;
		term2*=-zi*c;
		if( abs(c)<tol || abs(c)<tol*abs(sum)) break;
		if( abs(term)>abs(old))break;
		old=term;
		}
	return (sin(arg)*sum-cos(arg)*sum2)/sqrt(pi*sqrt(x));
/*printf(" Ai-=%e\n",(sin(arg)*sum-cos(arg)*sum2)/sqrt(pi*sqrt(x)));*/
	}
	}
if(z>0.)
	{
	x=sqrt(z);
	arg=.666666*x*z;
	return    .33333333*x*( imt(arg)-it(arg) );
/*return c1a*smallf(z)-c2a*smallg(z);*/
	}
x=sqrt(x);
arg=.666666*x*abs(z);
return .33333333*x*( jt(arg)+jmt(arg));
}


double bi(z) double(z);
{
int i,itmax=40;
double x,zeta,zi,pow(),sqrt(),arg,exp(),c,d,term,term2,old,old2,
	zz,tol=1.e-6,factor,sum,sin(),cos(),sum2;
if(z==0.) return c1a*sqrt(3.);
x=abs(z);
if(x>=asympt)
	{
	zeta= x*sqrt(x)*.6666666;
	zi=1./zeta;
	if(z>0.)
		{
		sum=0.;
		c=term=1.;
		old=1.;
		for(i=0;i<itmax;i++)
			{
			sum+=term;
			if(!i)factor=15.;
			else factor*= (6*i-5)*(6*i-3)*(6*i-1)/(2*i-1);
			c*= factor/(216.*(i+1));
			term*=zi*c;
			if( abs(c)<tol || abs(c)<tol*abs(sum)) break;
			if( abs(term)>abs(old))break;
			old=term;
			}
		return 1./sqrt(pi*sqrt(x))*exp(zeta)*sum;
/*printf(" Bi=%e\n", 1./sqrt(pi*sqrt(x))*exp(zeta)*sum);*/
		}
else
	{
	arg=zeta+pi*.25;
	sum=sum2=0.;
	d=c=term=1.;term2=-1.;
	old=old2=1.;
	for(i=0;i<itmax;i++)
		{
		if(i%2)sum2+=term2;
		else sum+=term;
		if(!i)factor=15.;
		else factor*= (6*i-5)*(6*i-3)*(6*i-1)/(2*i-1);
		c*=factor/(216.*(i+1));
		term*=-zi*c;
		term2*=-zi*c;
		if( abs(c)<tol || abs(c)<tol*abs(sum)) break;
		if( abs(term)>abs(old))break;
		old=term;
		}
	return (cos(arg)*sum+sin(arg)*sum2)/sqrt(pi*sqrt(x));
/*printf(" Bi-=%e\n", (cos(arg)*sum+sin(arg)*sum2)/sqrt(pi*sqrt(x)));*/
	}
	}
zz=sqrt(x);
arg=.666666*x*zz;
if(z>=0.)
return    zz*( imt(arg)+it(arg) )/sqrt(3.);
/*
return sqrt(3.)*(c1a*smallf(z)+c2a*smallg(z));*/
return zz/sqrt(3.)*(jmt(arg)-jt(arg));
}

double ibes(z,nu) double z,nu;
{
double pow(),exp(),sqrt(),gamma(),x,y,zz,sum,term,mu,ser,asym;
int k,n,itmax=500;
x=z*.5;
y=1.;
if(nu!=0.)y= pow(x,nu);
if(z<asympt)
	{
	sum=0.;
	zz=x*x;
	term=1.;
	for(k=0;k<itmax;k++)
		{
		asym=gamma(nu+1.+k);
		if(asym==0.)fprintf(stderr," pblm in ibes() gamma=%e nu=%e\n",asym,nu);
		sum+=term/asym;
		term*= zz/( ((double) (k+1)) );
		if( abs(term)<5.e-8 || abs(term)<1.e-8*abs(sum))return sum*y;
		}
fprintf(stderr," i bessel: precision not achieved, sum=%e\n",sum);
	return sum*y;
	}
/*else*/
mu=4.*nu*nu;
term=.125/z;/* for z>0 only:*/
return exp(z)/sqrt(2.*pi*z)*
		(1.+(mu-1.)*term*
		((.5-(mu-25.)*term*.1666666)*term*(mu-9.)-1.));
}

double dai(x) double x;
{
int k,itmax=40;
double zeta,zi,pow(),sqrt(),arg,exp(),c,d,term,term2,old,old2,
	tol=1.e-6,factor,sum,sin(),cos(),sum2;
double z,zz;
if(x==0.)return -c2a;
z=abs(x);
if(z>=asympt)
	{
	zeta= z*sqrt(z)*.6666666;
	zi=1./zeta;
	if(x>0.)
		{
		sum=0.;
		c=term=1.;
		old=1.;
		for(k=0;k<itmax;k++)
			{
			sum+=term;
			if(!k)factor=15.;
			else factor*= (6*k-5)*(6*k-3)*(6*k-1)/(2*k-1);
			c*=factor/(216.*(k+1));
			d=-(6*k+1.)/(6*k-1.)*c;
			term*=-zi*d;
			if( abs(c)<tol || abs(c)<tol*abs(sum)) break;
			if( abs(term)>abs(old))break;
			old=term;
			}
		return -.5*sqrt(sqrt(z)/pi)*exp(-zeta)*sum;
/*printf(" DAi=%e\n", -.5*sqrt(sqrt(z)/pi)*exp(-zeta)*sum);*/
		}
else
	{
	arg=zeta+pi*.25;
	sum=sum2=0.;
	d=c=term=1.;
	term2=-1.;
	old=old2=1.;
	for(k=0;k<itmax;k++)
		{
		if(k%2)sum2+=term2;
		else sum+=term;
		if(!k)factor=15.;
		else factor*= (6*k-5)*(6*k-3)*(6*k-1)/(2*k-1);
		c*=factor/(216.*(k+1));
		d=-(6*k+1.)/(6*k-1.)*c;
		term*=-zi*d;
		term2*=-zi*d;
		if( abs(c)<tol || abs(c)<tol*abs(sum)) break;
		if( abs(term)>abs(old))break;
		old=term;
		}
	return -(cos(arg)*sum+sin(arg)*sum2)*sqrt(sqrt(z)/pi);
/*printf(" Dai-=%e\n", -(cos(arg)*sum+sin(arg)*sum2)*sqrt(sqrt(z)/pi));*/
	}
	}
if(x<0.)
	{
	z=-x;
	zz=.666666*z*sqrt(z);
	return z*.333333*( jbes(zz,.666666)-jbes(zz,-.666666));
	}
zz=.666666*x*sqrt(x);
return -x*.333333*( ibes(zz,-.666666)-ibes(zz,.666666));
}

double dbi(x) double x;
{
double zeta,zi,pow(),sqrt(),arg,exp(),c,d,term,term2,old,old2,
	tol=1.e-6,factor,sum,sin(),cos(),sum2;
int k,itmax=40;
double z,zz;
if(x==0.)return c2a*sqrt(2.);
z=abs(x);
if(z>=asympt)
	{
	zeta= z*sqrt(z)*.6666666;
	zi=1./zeta;
	if(x>0.)
		{
		sum=0.;
		c=term=1.;
		old=1.;
		for(k=0;k<itmax;k++)
			{
			sum+=term;
			if(!k)factor=15.;
			else factor*= (6*k-5)*(6*k-3)*(6*k-1)/(2*k-1);
			c*=factor/(216.*(k+1));
			d=-(6*k+1.)/(6*k-1.)*c;
			term*= zi*d;
			if( abs(c)<tol || abs(c)<tol*abs(sum)) break;
			if( abs(term)>abs(old))break;
			old=term;
			}
		return sqrt(sqrt(z)/pi)*exp(zeta)*sum;
/*printf(" DBi=%e\n",sqrt(sqrt(z)/pi)*exp(zeta)*sum);*/
		}
else
	{
	arg=zeta+pi*.25;
	sum=sum2=0.;
	d=c=term=1.;term2=-1.;
	old=old2=1.;
	for(k=0;k<itmax;k++)
		{
		if(k%2)sum2+=term2;
		else sum+=term;
		if(!k)factor=15.;
		else factor*= (6*k-5)*(6*k-3)*(6*k-1)/(2*k-1);
		c*=factor/(216.*(k+1));
		d=-(6*k+1.)/(6*k-1.)*c;
		term*=-zi*d;
		term2*=-zi*d;
		if( abs(c)<tol || abs(c)<tol*abs(sum)) break;
		if( abs(term)>abs(old))break;
		old=term;
		}
	return (sin(arg)*sum-cos(arg)*sum2)*sqrt(sqrt(z)/pi);
/*printf("DBi-=%e\n", (sin(arg)*sum-cos(arg)*sum2)*sqrt(sqrt(z)/pi));*/
	}
	}
if(x<0.)
	{
	z=-x;
	zz=.666666*z*sqrt(z);
	return z*.57735027*(jbes(zz,-.666666)+jbes(zz,.666666));
	}
zz=.666666*x*sqrt(x);
return x*.57735027*(ibes(zz,-.666666)+ibes(zz,.666666));
}


double bigf(z) double z;
{
double pow(),sqrt(),x,y,zz,nu,sum,term,tk;
int k,n,itmax=5000;
sum=0.;
zz=z*z*z;
term=1.;
for(k=0;k<itmax;k++)
	{
	sum+=term;       tk=3.*k;
/*	term*= zz*((double)(3*k+1))/((double)((3*k+2)*(3*k+3)*(3*k+4)));*/
	term*= zz*(tk+1.)/((tk+2.)*(tk+3.)*(tk+4.));
	if( abs(term)<1.e-8 || abs(term/sum)<1.e-8) goto done;
	}
fprintf(stderr," bigf did not achieve desired precision\n");
done: return z*sum;
}
double bigg(z) double z;
{
double pow(),sqrt(),x,y,zz,nu,sum,term,tk;
int k,n,itmax=5000;
sum=0.;
zz=z*z*z;
term=.5;
for(k=0;k<itmax;k++)
	{
	sum+=term;tk=3.*k;
/*	term*= zz*((double)(3*k+2))/((double)((3*k+3)*(3*k+4)*(3*k+5)));  */
	term*= zz*((tk+2.))/(((tk+3.)*(tk+4.)*(tk+5.)));

	if( abs(term)<1.e-8 || abs(term/sum)<1.e-8)goto done;
	}
fprintf(stderr," bigg did not achieve desired precision\n");
done: return z*z*sum;
}

#define IAap  7.5
#define IAan 10.5
#define IBap 10.5
#define IBan 10.5

double IAi(z) double z;
{
double bigf(),bigg(),sqrt(),exp(),sin(),cos(),pow,x32,x34;
if( z==0.)return 0.;
if( z>0.)
	{
	if(z> IAap )
		{
		x32= z*sqrt(z);
		return 1./3.-.5/(sqrt(pi*x32))*exp(-.6666666*x32);
		}
	return c1a*bigf(z)-c2a*bigg(z);
	}
/*z<0.*/
if(z> - ( IAan ) )
	return -c1a*bigf(z)+c2a*bigg(z);
x32= -z*sqrt(-z);
return 2./3.-cos(.666666*x32+.25*pi)/sqrt(pi*x32);
}

double IBi(z) double z;
{
double sqrt(),exp(),sin(),cos(),pow,x32,x34;
if( z==0.)return 0.;
if( z>0.)
	{
	if(z>  IBap )
		{
		x32= z*sqrt(z);
		return exp(.6666666*x32)/(sqrt(pi*x32));
		}
	return sqrt(3.)*(c1a*bigf(z)+c2a*bigg(z));
	}
/*z<0.*/
if(z> - ( IBan ) )
	return -sqrt(3.)*(c1a*bigf(z)+c2a*bigg(z));
x32= -z*sqrt(-z);
return sin(.66666666*x32+.25*pi)/sqrt(pi*x32);
}

#define tolaw 1.e-8

aw(nu,z,jj,e) double z,nu,*jj,*e;
{
double x,zz,sum,term,term2,factor,
	arg,c,s,hnu,sum2,old,old2;
int k,l,n,itmax=500;

if(z>asympt)
	{
	x=1./z*z;	
	factor=nu*nu;
	sum=1.;
	sum2=1.;
	old=old2=1.;
	term=term2=1.;
	for(k=1;k<itmax;k+=2)
		{
		term*=-x*((k*k) -factor);
		if( abs(term)>abs(old))term=0.;
		old=term;
		sum+=term;
		l=k+1;
		term2*=-x*((l*l) -factor);
		if(abs(term2)>abs(old2))term2=0.;
		sum2+=term2;
		old2=term2;
		if( (abs(term)< tolaw || abs(term)< tolaw*abs(sum))
		&&(abs(term2)< tolaw || abs(term2)< tolaw*abs(sum2)))break;
		}
	sum2*=nu/z;
	*jj= jbes(z,nu)+ sin(nu*pi)/(pi*z)*(sum-sum2);
	zz=cos(nu*pi);
	*e=-ybessel-((1.+zz)*sum+(1.-zz)*sum2 )/(pi*z);
/* printf(" asympt %e %e\n",*jj,*e);*/
	 return 0;
	}
x=z*.5;
hnu=.5*nu;
sum=sum2=0.;
zz=x*x;
factor=1.;
for(k=0;k<itmax;k++)
	{term=factor/(gamma(1.+k-hnu)*gamma(hnu+1.+k));
	sum+=term;
	term2=factor/(gamma(1.5+k-hnu)*gamma(hnu+1.5+k));
	sum2+=term2;
	factor*= -zz;
	if( abs(term2)<tolaw*abs(sum2))goto d1;
	}
fprintf(stderr," tolerance not met aw series %le %le\n",term2,sum2);
d1:sum2*=x;
arg= pi*hnu;
c=cos(arg);s=sin(arg);
*jj= c*sum+s*sum2;
*e= s*sum-c*sum2;
return 0;
}
