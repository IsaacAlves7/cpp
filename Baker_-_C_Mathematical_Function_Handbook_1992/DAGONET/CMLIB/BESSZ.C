/*  Zeros of bessel functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

zerobes(a,n,z,d) double a,z[];int d,n;
	where a is order, z[] are first n zeros,
	d= 1 for J 2 for Y 3 for J' 4 for Y'
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"


#define pio2 pi*.5

double fi(y) double y;
{/* solve tan(fi)-fi=y for fi*/
double p,pp,q,r,pow(),atan();
if(y==0.)return 0.;
if(y>1.e5)return pio2;
if(y<1.)
	{
	p=pow(3.*y,.3333333333);
	pp=p*p;
	p*=(1.+pp*(-210.+pp*(27.-2.*pp)))/1575.;
	}
else
	{
	p=1./(y+pio2); pp=p*p;
	p=pio2-p*(1.+pp*(2310.+pp*(3003.+pp*(4818.+pp*(8591.+pp*16328)))))/3465.;
	}
q=y+p;
pp=q*q;
r=(p-atan(q))/pp;
return p-(1.+pp)*r*(1.+r/q);
}


double ybessl,sum1,sum2;

double jass(z,nu) double z,nu;
{
double term2,term,x,y,q,arg,c,s;
double zz,f1,f2,f1o,f2o,nuph,numh,npho,mu;
int k,twok,n,itmax=30;
/* INCREASING NU ONLY MAKES IT WORSE*/
/*if( abs(nu)<1.) return jass(z,nu+1.)*2.*(nu+1.)/z-jass(z,nu+2.);*/

z=abs(z);

y= sqrt(2./(pi*z));
sum1=0.;
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
		sum1+=f1;
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
	if( (abs(term)<1.e-7 || abs(term/sum1)<1.e-6)&&
	 ( abs(term2)<1.e-7 || abs(term2/sum2)<1.e-6)  )break;
	}
arg= z-pi*(nu*.5+.25);
return y*(cos(arg)*sum1-sin(arg)*sum2);
*/
/*HANKEL EXPANSION:*/
zz=1./(64.*z*z);
term=(mu-1.)*(mu-9.)*zz;
sum1=1.-term*.5+term*(mu-25.)*(mu-49.)*zz/24.;
term2=(mu-1.)*.125/z;
sum2=term2*(1.-(mu-9.)*(mu-25.)*zz/6. );
arg= z-pi*(nu*.5+.25);
c=cos(arg);s=sin(arg);
ybessl=y*(s*sum1+c*sum2);
return y*(c*sum1-s*sum2);
}

besspq(a,x,pa,qa,pa1,qa1) double a,x, *pa,*qa,*pa1,*qa1;
{
double jv,yv,s,c,chi,garb,ck;
struct complex z,jc,yc,h2,jp,yp,hp;
int n,ivk;
n=a;            z.x=x;z.y=0.;
if(x>10.)
	{
	jass(x,a); *pa=sum1;*qa=sum2;
	jass(x,a+1.); *pa1=sum1;*qa1=sum2;
	}
else
	{
	if( (a- (double)n)==0.)
		{
		bessel( n, &z,&jc,&yc,&h2,&jp,&yp,&hp,&ivk);
		jv=jc.x;yv=yc.x;
		chi=x-pi*(.25+.5*a);
		s=sin(chi);c=cos(chi);
		garb=sqrt((pi*x)*.5);
		*pa= garb*(c*jv+s*yv);
		*qa= garb*(c*yv-s*jv);
		ck= (*pa*c-*qa*s);
		n++;
		bessel( n, &z,&jc,&yc,&h2,&jp,&yp,&hp,&ivk);
		jv=jc.x;yv=yc.x;
		chi=x-pi*(.75+.5*a);
		s=sin(chi);c=cos(chi);
		*pa1= garb*(c*jv+s*yv);
		*qa1= garb*(c*yv-s*jv);
		}
	else
		{/* order not integer, not in asymptotic regime for jas*/
		jv=jbes(x,a);yv=ybes(x,a);
		chi=x-pi*(.25+.5*a);
		s=sin(chi);c=cos(chi);
		garb=sqrt((pi*x)*.5);
		*pa= garb*(c*jv+s*yv);
		*qa= garb*(c*yv-s*jv);
		ck= (*pa*c-*qa*s);
		jv=jbes(x,a+1.);yv=ybes(x,a+1.);
		chi=x-pi*(.75+.5*a);
		s=sin(chi);c=cos(chi);
		*pa1= garb*(c*jv+s*yv);
		*qa1= garb*(c*yv-s*jv);
		}
	}
return 0;
}

zerobes(a,n,z,d) double a,z[];int d,n;
{
double e=1.e-4;/* desired accuracy. should be consistent with
accuracy of jas p,q values*/
double ck,ckd;
double aa,a1,a2,b,bb,c,chi,co,mu,mu2,mu3,mu4,p,pa,pa1,pp1,p0,p1,q1,
psi,q,qa,x,qq1,si,t,tt,u,v,w,xx,x4,y,sqrt(),cos(),pow(),r0,qa1;
int j,s;
aa=a*a; mu=4.*aa;mu2=mu*mu;mu3=mu*mu2;mu4=mu2*mu2;
if(d<3)
	{
	p=7.*mu-31.;p0=mu-1.;
	if(1.+p == p) p1=q1=0.;
	else
		{
		p1=4.*(253.*mu2-3722.*mu+17869.)*p0/(p*15.);
		q1=1.6*(83.*mu2-982.*mu+3779.)/p;
		}
	}
else
	{
	p=7.*mu2+82.*mu-9.; p0=mu+3.;
	if(1.+p == p) p1=q1=0.;
	else
		{
		p1=(4048.*mu4+131264.*mu3-221984.*mu2-417600.*mu+1012176.)/(60.*p);
		q1=1.6*(83.*mu3+2075.*mu2-3039.*mu+3537.)/p;
		}
	}
t= (d==1 || d==4)? .25 : .75;
tt=4.*t;
if(d<3)
	{pp1=5./48.;qq1=-5./36.;}
else
	{pp1=-7./48.;qq1=35/288.;}
y=.375*pi;bb= (a>=3.)? pow(a,-.6666666): 1.;
a1=3*a-8.;psi=pi*(.5*a+.25);
for(s=1;s<=n;s++)
	{
	if( s==1 && a==0. && d==3)
		{x=0.;j=0;}
	else
		{
		if( s>=a1 )
			{
			b=(s+.5*a-t)*pi;c=.0125625/(b*b);
			x=b-.125*(p0-p1*c)/(b*(1.-q1*c));
			}
		else
			{
			if(s==1)
				{
				switch (d)
					{
					case 1 : x=-2.33811;break;
					case 2 : x=-1.17371;break;
					case 3:  x=-1.01879;break;
					default: x=-2.29444;break;
					}
				}
			else
				{
				x=y*(4.*(s)-tt);v=1./(x*x);
				x= -pow(x,.6666666)*(1.+v*(pp1+qq1*v));
				}
		u=x*bb;
		v= fi(.666666*pow((-u),1.5));
		w=1./cos(v);
		xx=1.-w*w;
		c=sqrt(u/xx);
		x=w*(a+c/(48.*a*u)) *((d<3)? -5./u-c*(-10./xx+6.):
		( 7./u+c*(-14./xx+18.)));
		}
		j=0;
		do
			{
			xx=x*x; x4=xx*xx; a2=aa-xx;
		besspq(a,x,&pa,&qa,&pa1,&qa1);
			chi=x-psi;
			si=sin(chi);co=cos(chi);
			switch (d)
				{
				case 1: r0= (pa*co-qa*si)/(pa1*si+qa1*co);break;
				case 2: r0= (pa*si+qa*co)/(qa1*si-pa1*co);break;
				case 3:	
					r0=a/x-(pa1*si+qa1*co)/(pa*co-qa*si);break;
				default:
					r0=a/x-(qa1*si-pa1*co)/(pa*si+qa*co);break;
				}
			j++;
			if(d<3)
				{
				u=r0;w=6.*x*(2.*a+1.);
				p=(1.-4.*a2)/w;
				q=(4*(xx-mu)-2.-12.*a)/w;
				}
			else
				{
				u=-xx*r0/a2;
				v=2.*x*a2/(3.*(aa+xx));
				w=64.*a2*a2*a2;
				q=2.*v*(1.+mu2+32.*mu*xx+48.*x4)/w;
				p=v*(1.+(-mu2+40.*mu*xx+48.*x4)/w);
				}
			w=u*(1.+p*r0)/(1.+r0*q);x+=w;
			} while( abs(w/x)>e && j<5);
		}/*else of a=0,d==3,s==1 */
	z[s-1]=x;
	}
return 0;
}
