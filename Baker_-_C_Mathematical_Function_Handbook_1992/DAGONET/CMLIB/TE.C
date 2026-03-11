/* basic elliptic functions real arguments
	copyright 1991 Louis Baker. All rights reserved.

	tek()	returns complete elliptic integrals 1st,2nd kind
		id=0 for 2nd argument m else =1 for m1 supplied
		0<=m<=1

	tef()	incomplete elliptic integrals of 1st,2nd kind
		as functions of phi,m. see e3m.c for third kind
		any real m

	jzeta	Jacobian Capital Zeta function

	am() 	amplitude

	jef()	Jacobian elliptic function sn,cn,dn all real m
		returns in global double amplitude am

	heuman	Heuman's Lambda function

	theta	Jacobian theta functions i=1,4 real u,m

	neville	Neville theta functions

tek,tef,jef jzeta based upon routines in AFWL math notes.
Modified to handle  m<0 and m>1, m==0,1 as special cases,
calculate Capital Theta, return amplitude from jef
*/


#include "cmlib.h"
#include "protom.h"
static teff(double phi,double m,double sig,double *f,double *e);

double kimag,eimag;

tek(id,m,k,e)
double m,*k,*e;int id;
{/* complete elliptic integrals of 1st,2nd kinds*/
int i,n;
double sqrt(),log(),
	rpk,rkn,rkp[60],tol=1.e-6,pk2,pkp,gol,gk,fk,fe,ge,c,d,h,t1,t2,aux;
kimag=eimag=0.;
if(id) rpk=sqrt(m);
if(!id)aux=m;else aux=1.-m;
if(aux>1. )
	{
	aux=1./aux;
	h=sqrt(aux);d=(1.-aux);
	tek(0,aux, &fk,&fe);
	tek(0,d, &gk,&ge);
	*k=h*fk;kimag=h*gk;
	*e=(fe-fk*d)/h;eimag=(ge-aux*gk)/h;		
	return 0;
	}
if(m==1.)
	{*k=1.e60; *e=1.;return 0;}
*k=pi*.5;
*e=*k;
if(aux<0.)
	{
	aux= -aux/(1.-aux);
	h=sqrt(1.-aux);
	tek(0,aux,k,e);
	(*k) *=h; (*e)/=h;
	return 0;
	}
if(m==0.)return 0;
if(m> .999)
	rpk=sqrt(1.-m);
else
	{
	rkn=sqrt(m);
	for(i=0;i<60;i++)
		{
		rkp[i]=aux= sqrt(1.-rkn*rkn);
		rkn= (1.-aux)/(1.+aux);
		if(i && rkn<tol)break;
		}
	n=i;
	for(i=n;i>=0;i--)
		{
		rkn=rkp[i];
		aux=1.+rkn;
		*k *= 2./aux;
		*e= aux* *e -*k *rkn;
		}
	return 0;
	}
/* m>.999 */
pk2=rpk*rpk;
pkp=pk2;
gol= log(4./rpk);
gk=gol-1.;
fk=.25;
fe=.25;
*k=gol+fk*gk*pkp;
*e=1.+.5*(gol-.5)*pkp;
ge=gk;
for(i=2;i<2000;i++)
	{
	aux=i+1;
	d=i;
	pkp=pkp*pk2;
	c=d/aux;
	*k *= d*d/(aux*aux);
	*e *=c;
	h=1./(d*aux);
	gk -= 1./(d*d);
	ge-=h;
	t1= *k *gk*pkp;
	*k += t1;
	t2= fe*ge*pkp;
	*e += t2;
	if(t1<tol && t2<tol)return 0;
	fe *=c;
	ge -=h;
	}
return 0;
}

tef(phi,m,sig,f,e) double phi,m,sig,*f,*e;
{ double sin(),log(),tan(),sqrt(),sin(),asin(),am();
double  u,ff,ec,w,plus,plus1,t,u1,a,p,r,bd,ek,ee;
struct complex carg,cans;
if(m==0.)
	{*e=*f=phi;return 0;}
if(m==1.)
	{
	*e=sin(phi);*f= log(tan(pi*.25+.5*phi));
	return 0;
	}
if(m>1.)
	{/* A&S 17.4.15-16*/
	plus= sqrt(m);plus1=1./m;
	ff=plus*sin(phi);
	if(ff>1. || ff<-1.)
		{/* imaginary angle theta- need complex functions*/
		*f=*e=errorcode;return 1;
		}
	t= asin(ff);
	u=sig;
	tef(t,plus1,u,&ff,&ec);
	*f= ff/plus; /* *f= u */
	w=*f*plus; p=am(w,plus1);/*am plus1 or m???*/
	tef(p,plus1,sig,&ff,&ec);
	*e= plus*ec-(m-1.)* *f;
	return 0;
	}
else if(m<0.)
	{/* A&S 17.4.17 F, homebrew E note A&S use -m for m */
	plus=1./(1.-m);plus1=-m*plus;
	tek(0,plus1,&ek,&ee);
	tef(.5*pi-phi,plus1,sig,&ff,&ec);
	r=sqrt(plus);bd=1./r;
	*f= r*( ek-ff);
	w= *f *bd;
	u=sqrt(plus1);u1=sqrt(1.-plus1);
	r=sin(phi);a= u*r;ek=sqrt(u1*u1+ a*a);
	a= r/ek;
	bd=asin(a);
	tef(bd,plus1,sig,&ff,&ec);p=u*a;
	*e= (ec-plus1*a*cos(bd)/sqrt(1.-p*p))/u1;
	return 0;
	}
teff(phi,m,sig,f,e);
return 0;
}

static teff(phi,m,sig,f,e) double phi,m,sig,*f,*e;
{/* CAVEAT- USER SHOULD CALL TEF, NOT TEFF*/
double sqrt(),tan(),atan(),sin(),cos(),asin(),log(),tol=1.e-8;
double phii,ek,ee,plus,plus1,w,signem,ph,rk,ss,sk,alphar;
double d,g,ps,u1,u,h,a,ap,cnk,be,reler,oreler,p,or,pr,bf,s,c;
double sumem,bd,r,t,p2,t1,t2,pk,ce,ff,am(),b,lou;
double aa[50],bb[50],cc[50],psav[50];
int i,j,n,nk,stop,ios,ins,i4,nq,m2p,it,k;
if(phi<0.)
	{
	w=-1.;
	ph=-phi;
	}
else
	{
	w=1.;
	ph=phi;
	}
rk=sqrt(m);
n=ph/(2.*pi);
a=ph- n*2.*pi;
b=a/pi *2.;
k=b;
nq=k+1;
switch (nq)
	{
	case 1: nk= (n<<2);
		signem=1.;
		ap=a;
		break;
	case 2:
		nk= (n<<2)+2;
		signem=-1.;
		ap=pi-a;
		break;
	case 3:
		nk=(n<<2)+2;
		signem=1.;
		ap=a-pi;
		break;
	default:
		nk= (n<<2)+4;
		signem=-1;
		ap= 2.*pi-a;
		break;
	}
cnk=nk;
phii=ap;
tek(0,m,&ek,&ee);
plus=cnk*ek;
plus1=cnk*ee;
it=0;
if(abs(phii-.5*pi)<tol)it=1;
if(abs(rk-1.)<tol)
	{
	it++;
	if(it==2)
		{
		*f=w*1.e70; *e= w*(plus1+signem);return 0;
		}
	*f=w*(plus+signem*log(tan(phii*.5+.785398163397448)));
	*e=w*(plus1+signem*sin(phii));
	return 0;
	}
if(abs(rk)<tol)
	{
	*f= w*(plus+signem*phii);
	*e=w*(plus1+signem*phii);
	return 0;
	}
it++;
if(it==2)
	{
	tek(0,m,&ek,&ee);
	*f=w*(plus+signem*ek);
	*f=w*(plus1+signem*ee);
	return 0;
	}
if(abs(phii)<1.e-30)
	{
	*f=w*plus;
	*e=w*plus1;
	return 0;
	}
if(m<=.75)
{
tek(0,m,&ek,&ee);
s=sin(phii);
c=cos(phii);
sk=m;
ce=2.*phii/pi;
t2=ce*ek;
t1=ce*ee;
a=.5;
t=.5*a*sk;
r=t;
ss=s*s;
ps=1.;
h=.5;
ff=.5;
pk=sk;
u1=10.;
for(i=2;i<20000;i++)
	{
	j=i<<1;
	d= j-1;
	g=j-3;
	ee=1./j;
	ps *=ss;
	a=ee*(d*a+ps);
	ff*=d*ee;
	h*=g*ee;
	pk*=sk;
	u=ff*a*pk;
	if( u1*u1/(u1-u) < sig)break;
	u1=u;
	t+=u;
	r+=h*a*pk;
	}
(*f)=w*((t2-s*c*t)*signem+plus);
(*e)=w*((t1+s*c*r)*signem+plus1);
return 0;
}
	/* else m>=.75   */
	alphar=asin(rk);
	aa[0]=1.;
	bb[0]=cos(alphar);stop=49;
	for(i=1;i<50;i++)
		{j=i-1;
		aa[i]=.5*(aa[j]+bb[j]);
		bb[i]=sqrt(aa[j]*bb[j]);
		cc[i]=.5*(aa[j]-bb[j]);
		if(abs(cc[i])<sig){stop=i;break;}
		}
	p=phii;
	p2=1.;
	nq=0;
	ios=1;
	m2p=0;
	i4=0;
	oreler=1.e20;
	or=oreler;
	for(i=0;i<=stop;i++)
		{
		psav[i]=p;
		p2 *=2.;
		bd= tan(p)*bb[i]/aa[i];
		bf=atan(bd);
		label1:
			ins= (bf<0.) ? -1 : 1 ;
		if( ios*ins < 0)
			{
			nq++;
			nq%=4;
    		}
		switch (nq)
			{
			case 0:
				if(i4)
					{
					i4=0;
					m2p++;
					}
				be=bf+(2*m2p)*pi;
				break;
			case 1:
			case 2:
				be=bf+(2*m2p)*pi+pi;
				break;
			case 3:
				be=bf+(2*m2p)*pi+2.*pi;
				i4=1;
			}
		ios=ins;
		pr=p/be;
		reler= abs(or-pr)/(pr+or);
		if( oreler<reler )
				{
				ios=-ios;
				goto label1;
				}
		p +=be;
		or=pr;
		oreler=reler;
		}
	lou=(p/(p2*aa[stop]));
	*f= w*(plus+signem*lou);
	tek(0,m,&ek,&ee);
/*printf(" debug K,E=%le %le lou=%le\n",ek,ee,lou);*/
	sumem=0.;
	for(i=1;i<=stop;i++)
		sumem += cc[i]*sin(psav[i]);
/* not F(phi) but transformed phi needed here? */
	*e=w*(plus1+signem*(ee/ek* (lou) + sumem));
/*	INCORRECT CODE FROM MN70:
	*e=w*(plus1+signem*(ee/ek* *f + sumem));*/
/*printf(" E debug: plus1=%le sumem=%le F=%le, sigmen=%le\n"
,plus1,sumem,*f,signem);*/
	return 0;
}

double jzeta(phi,m) double phi,m;
{/*jacobian zeta function*/
double k,e,f,ee,sig=1.e-5;
tek(0,m,&k,&e);
tef(phi,m,sig,&f,&ee);
if(ee==errorcode || f==errorcode)return errorcode;
return ee-f*e/k;
}

double am(u,m) double u,m;
{/* returns phi= am u inverse function of u=F(phi,m=sin^2 alpha) */
/* asin() returns within range -pi/2 to pi/2  am(u+2K)= pi+am(u) */
double sn,cn,dn,asin(),mult,e,k,offset;
offset=0.;mult=1.;
if(u<0.){u=-u;mult=-1.;}
tek(0,m,&k,&e);e=((int)(u/(2.*k)));
u-=e;offset=pi*e;
jef(u,m,&sn,&cn,&dn);
return  mult*(asin(sn)+offset);
}


double jtheta,amplitude;

jef(u,m,sn,cn,dn) double u,m,*sn,*cn,*dn;
{
int i,n;
double a[200],c[200],phi[200],v,am,am1,b,argu,t,
	asin(),sin(),cos(),sqrt(),pow(),tol=1.e-7;

/* Jacobian theta function valid for all real m<=1 */
double exp(),log(),twon,term,k,sum,lou,sqrtm1,mu,mu1;

if(m<0.)
	{/* A&S 16.10.1-4*/
	t=-m;mu= t/(1.+t);mu1=1./(1.+t);k=sqrt(mu1);v=u/k;
	jef(v,mu,&term,&sum,&lou);
	*dn= 1./lou;
	*sn= k* term* *dn;
	*cn= sum* *dn;
	return 0;
	}
if(m>1.)
	{ /* A&S 16.11.1-4*/
	mu=1./m; k=sqrt(mu);v=u/k;
	jef(v,mu,&term,&sum,&lou);
	*sn=k*term;
	*cn=lou;
	*dn=sum;
	return 0;
	}
v=u;am=m;
if(m==0.)
	{
	*sn=sin(u);
	*cn=cos(u);
	*dn=1.;
	return 0;
	}
else if(m==1.)
	{
	*sn=tanh(u);
	*cn=1./sinh(u);
	*dn=*cn;
	return 0;
	}
am1=1.-am;
a[0]=1.;
sqrtm1=sqrt(am1);
b=sqrtm1;
/*c[0]=sqrt(m);not used anywhere*/
for(i=1;i<200;i++)
	{
	a[i]=.5*(a[i-1]+b);
	c[i]=.5*(a[i-1]-b);
	if(abs(c[i])<tol)break;
	b=sqrt(b*a[i-1]);
	}
n=i;
twon=pow(2., (double)(n));
phi[n]=a[n]*v* twon;

sum=0.;
term= .5/twon;
for(i=n;i>0;i--)
	{
	argu=c[i]*sin(phi[i])/a[i];
	t=asin(argu);
	phi[i-1]=.5*(t+phi[i]);
	sum-= term*log(cos(2.*phi[i-1]-phi[i]));
	term*=2.;
	}
argu=phi[0]; amplitude=argu;
*sn=sin(argu);
*cn=cos(argu);
lou= cos(phi[1]-argu);
*dn= *cn/lou;
tek(0,m,&k,&term);
jtheta=sum+.5*log(2.*sqrtm1/pi*k * lou/cos(argu));
jtheta=exp(jtheta);
return 0;
}

double heuman(phi,m) double phi,m;
{
double k,e,ff,ee;
tek(0,m,&k,&e);
tef(phi,1.-m,1.e-3,&ff,&ee);
return 2./pi*(k*ee-(k-e)*ff);
}

theta(v,m,t1,t2,t3,t4)
double v,m,*t1,*t2,*t3,*t4;
{
double sn,cn,dn,bigk,k,kp,e,u,srk,srkp;
tek(0,m,&bigk,&e);
u=2.*bigk*v/pi;
jef(u,m,&sn,&cn,&dn);
k=sqrt(m);kp=sqrt(1.-m);srk=sqrt(k);srkp=sqrt(kp);
*t1= sn*srk*jtheta;
*t2=cn*srk/srkp*jtheta;
*t3=dn/srkp*jtheta;
*t4=jtheta;
return 0;
}

neville(u,m,ts,tc,td,tn) double u,m,*ts,*tc,*tn,*td;
{
double v,e,bigk,t1,t2,t3,t4,t10,t20,t30,t40;
tek(0,m,&bigk,&e);
v=u*pi/(2.*bigk);
theta(v,m,&t1,&t2,&t3,&t4);
theta(0.,m,&t10,&t20,&t30,&t40);
*ts=t1*bigk*2./(pi*t20*t30*t40);
*tc=t2/t20;
*td=t3/t30;
*tn=t4/t40;
return 0;
}