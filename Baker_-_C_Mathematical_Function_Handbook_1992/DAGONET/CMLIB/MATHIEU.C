/*  Mathieu functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

based upon FORTRAN code of CALGO 352 of D. S. Clemm

*/

#include <stdio.h>
#include <alloc.h>
#include "cmlib.h"
#include "protom.h"

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))
#define INDEX(i,j) [(j)+(i)*(coln)]

static double qmath,tol,qinv;
static double a,a0,a1;
static int type;

static double gmath[2][200],dg[2][200];
static int mf,m0,m1,m2s;
#define fl 1.e30

/* the following prototypes are not intended to be
global, as they have fairly common names.
comment the prototypes out if your C compiler is not
ANSI compatible
*/
static int coef( void);

int tmofa(alfa,tm,dtm) double alfa,*tm,*dtm;
{int k,kk,kt,l;
/* NOT modified to do c.f. the Lenz/Thompson/Barnett method*/
double aa,dtype,
q1,q2,t,tt,sqrt(),a[3],b[3],aux;

kt=0;aa=alfa;dtype=type;
for(l=0;l<2;l++) {for(k=0;k<200;k++)
			{gmath[l][k]=dg[l][k]=0.;}
		}
if( (type%2))m0=3;
else m0=type+2;
k=.5*sqrt(max(0.,3.*qmath+aa));
m2s=min(2*k+m0+4,398+m0%2);
a[0]=1.;a[1]= (aa-(m2s+2.)*(m2s+2.))*qinv;
b[0]= (aa-m2s*m2s)*qinv;
b[1]=a[1]*b[0]-1.;
q1=a[1]/b[1];
kt=1;
for(k=0;k<200;k++)
	{
	mf=m2s+4+2*k;
	t=(aa-mf*mf)*qinv;
	a[2]=t*a[1]-a[0];
	b[2]=t*b[1]-b[0];
	q2=a[2]/b[2];
	if(abs(q1-q2)<tol){kt=0;break;}
	q1=q2;
	a[0]=a[1];b[0]=b[1];
	a[1]=a[2];b[1]=b[2];
	/* no need to renormalize???*/
	}
t=1./t;
tt=-t*t*qinv;
l=mf-m2s;
for(k=2;k<=l;k+=2)
	{
	t=1./(qinv*(aa-(mf-k)*(mf-k))-t);
	tt=t*t*(tt-qinv);
	}
kk= m2s>>1;
if(kt==1)q2=t;
gmath[1][kk]=.5*(q2+t);
dg[1][kk]=tt;
gmath[0][1]=1.;
for(k=m0;k<=m2s;k+=2)
	{
	kk=k>>1;
	if(k<5)
		{
		if(k<3)
			{
			gmath[0][1]=aa*qinv;
			dg[0][1]=qinv;
			}
		else if(k==3)
			{
			gmath[0][1]=(aa-1.)*qinv+dtype-2.;
			dg[0][1]=qinv;
			}
		else
			{
			aux=gmath[0][1];
			gmath[0][2]=(aa-4.)*qinv+(dtype-2.)/aux;
			dg[0][2]=qinv+(2.-dtype)*dg[0][1]/(aux*aux);
			if(type==2)gmath[0][1]=0.;
			}
		}
	else
		{
		aux=gmath[0][kk-1];
		gmath[0][kk]=(aa-(k-2)*(k-2))*qinv-1./aux;
		dg[0][kk]=qinv+ dg[0][kk-1]/(aux*aux);
		}
	if(abs(gmath[0][kk])<1.)goto stage2;
	}
/* backtrack*/
*tm= gmath[1][kk]-gmath[0][kk];
*dtm=dg[1][kk]-dg[0][kk];
m1=m2s;
kt=m2s-m0;
for(l=2;l<=kt;l+=2)
	{
	k=m2s-l;
	kk=k>>1;
	gmath[1][kk]=1./(qinv*(aa-k*k)-gmath[1][kk+1]);aux= gmath[1][kk];
	dg[1][kk]=-aux*aux*(qinv-dg[1][kk+1]);
	if(k<=2)
		{
		gmath[1][1]*=2.;
		dg[1][1]*=2.;
		}
	t= gmath[1][kk]-gmath[0][kk];
	if(abs(t)<abs(*tm))
		{*tm=t;*dtm= dg[1][kk]-dg[0][kk];m1=k;}
	}
stage2:
m1=k; k=m2s;kk=k>>1;
redo:
	if(k==m1)
	{
	if(k<=2)
		{gmath[1][1]*=2.;dg[1][1]*=2.;}
	*tm=gmath[1][kk]-gmath[0][kk];
	*dtm=dg[1][kk]-dg[0][kk];
	return 0;
	}
k-=2;
kk--;
t= (aa-k*k)*qinv-gmath[1][kk+1];
if(abs(t)>=1.)
	{again:
	gmath[1][kk]=1./t;
	dg[1][kk]=( dg[1][kk+1]-qinv)/(t*t);
	goto redo;
	}
if(k==m1)
	{
	if(t==0.)return 1;
	goto again;
	}
q1= dg[1][kk+1]-qinv;
while(1)
{
gmath[1][kk]=fl;
gmath[0][kk]=t;
k-=2;
kk--;
q2= t*(aa-k*k)*qinv-1.;
if(k==m1)
	{
	if(q2==0.)return 1;
	gmath[1][kk]=t/q2;
	dg[1][kk]= (q1-qinv*t*t)/(q2*q2);
	goto redo;
	}
if(abs(q2)-abs(t)<0.)
	{q1=q1/(t*t)-qinv;
	t=q2/t;
	continue;
	}
else/* was goto 280*/ 
	{
	gmath[1][kk]=t/q2;
	dg[1][kk]= (q1-qinv*t*t)/(q2*q2);
	goto redo;
	}
}/*stage3 loop*/
}

static int l,s,p,klast,kmax,mm,ml,n;
static double a,dmax,dlast,total,t,ab[200],u1,u2,x;

static int coef()
{int k,ka,kb,kk,m;
double log(),v2=1.e-15;
m=tmofa(a,&t,&t);/*set mf,etc and g*/
if(m){fprintf(stderr," failure in tmfoa from coef\n");return m;}
for(k=0;k<200;k++)ab[k]=0.;
ka=m1-m0+2;
for(k=2;k<=ka;k+=2)
	{
	kk=(m1-k)>>1;
	if(k<=2)ab[kk]=1.;
	else ab[kk]=ab[kk+1]/gmath[0][kk+1];
	}
ka=0;
for (k=m1;k<=m2s;k+=2)
	{
	kk=k>>1;
	ml=k;
	if( gmath[1][kk]!=fl)ab[kk]=ab[kk-1]*gmath[1][kk];
	else
		{
		t=ab[kk-2];
		if(k==4 && m1==2)t+=t;
		ab[kk]=t/( (a-(k-2)*(k-2))*qinv*gmath[0][kk]-1.);
		}
	if( abs(ab[kk])>=1.e-12)ka=0;
	if(ka==5)goto norm;
	ka++;
	}
t=log(abs(ab[kk])/v2)/ log(1./abs(gmath[1][kk]));
ka=2*t;
ml=ka+2+m2s;
if(ml>399)
	{
	fprintf(stderr," ml ka m2s %d %d %d %e\n",ml,ka,m2s,t);
	return -1;
	}
kb=ka+2+mf;
t=1./(qinv*(a-kb*kb));
kk=mf-m2s;
for(k=2;k<=kk;k+=2) t=1./(qinv*(a-(kb-k)*(kb-k))-t);
kk=(ml>>1);
gmath[1][kk]=t;
for(k=2;k<=ka;k+=2)
	{kk=(ml-k)>>1;
	gmath[1][kk]=1./( qinv*(a-(ml-k)*(ml-k))-gmath[1][kk+1]);
	}
ka=m2s+2;
for(k=ka;k<=ml;k+=2)
	{
	kk=k>>1;
	ab[kk]=ab[kk-1]*gmath[1][kk];
	}
norm: t=ab[0];
mm= type%2; ka=mm+2;
for(k=ka;k<=ml;k+=2)
	{kk=k>>1;
	if( abs(t)<abs(ab[kk]))
		{
		t=ab[kk];mm=k;
		}
	}
for(k=0;k<=kk;k++) ab[k]/=t;
return 0;
}

/* fj returns J(n)  for z= u1; fy returns J or Y  z=u2, depending on sol
dy,dj return the corresponding derivatives*/
double ju1[200],jyu2[200],dju1[200],djyu2[200];


double fj(n) int n;
{return ju1[n];
}

double fy(n) int n;
{return jyu2[n];
}

double dy(n) int n;
{return djyu2[n];
}

double dj(n) int n;
{return dju1[n];
}



double ds(kk)int kk;
{
double dsv;
int k,n1,n2;
k=kk;
n1=k-s;n2=k+s+p;
dsv= ab[k]*( fj(n1)*fy(n2)-fj(n2)*fy(n1));

/*printf(" ds %d %d %d %le %le %le %le %le\n"
	,k,n1,n2,ab[k],fj(n1),fy(n2),fj(n2),fy(n1));
*/
if( (k+n)%2)dsv=-dsv;return dsv;
}

double dc(kk)int kk;
{
double dsv;
int k,n1,n2;
k=kk;
n1=k-s;n2=k+s+p;
dsv= ab[k]*( fj(n1)*fy(n2)+fj(n2)*fy(n1));
/*printf(" dc %d %d %d %le %le %le %le %le\n"
	,k,n1,n2,ab[k],fj(n1),fy(n2),fj(n2),fy(n1));
*/
if((s+p)==0.)dsv*=.5;
if( (k+n)%2)dsv=-dsv;return dsv;
}


double dds(kk)int kk;
{
double dsv;
int k,n1,n2;
k=kk;
n1=k-s;n2=k+s+p;
dsv= ab[k]*(u2*fj(n1)*dy(n2)-fj(n2)*dy(n1))-u1*(fy(n2)*dj(n1)-fy(n1)*dj(n2));
if( (k+n)%2)dsv=-dsv;return dsv;
}
double ddc(kk)int kk;
{
double dsv;
int k,n1,n2;
k=kk;
n1=k-s;n2=k+s+p;
dsv= ab[k]*(u2*fj(n1)*dy(n2)+fj(n2)*dy(n1))-u1*(fy(n2)*dj(n1)+fy(n1)*dj(n2));
if((s+p)==0.)dsv*=.5;
if( (k+n)%2)dsv=-dsv;return dsv;
}

double ps(k)
{double sin();
return ab[k]*sin( x*(p+(k<<1)));
}

double pc(k)
{double cos();
return ab[k]*cos( x*(p+(k<<1)));
}

double dps(k)
{double cos();
t=p+(k<<1);
return ab[k]*t*cos( x*t);
}
double dpc(k)
{double sin();
t=p+(k<<1);
return -ab[k]*t*sin( x*t);
}


#define normlim  1.e-14

sum(func) double (*func)();
{       int k;
k=0; total= func(0); dmax=total;t=abs(total); kmax=0;
if(t<normlim)/* for x= pi/2 when f should be zero*/
	{
	kmax=klast=0;
	dmax=total=0.;
	t=1.; /* prevent division by zero*/
	return 1;
	}
/*printf(" in sum l=%d\n",l);*/
for(klast=1;klast<=l;klast++)
	{
	dlast= func(klast);
	/*printf(" func(%d=klast)=%le\n",klast,dlast);*/
	total+=dlast;
	if(t<abs(dlast))
		{
		dmax=dlast;
		t=abs(dmax);
		kmax=klast;
		}
	if(klast>s)
		{/* stop if "many" (3) small terms*/
		if( abs(dlast)/t> 1.e-12)k=0;
		k++;
		if(k>=3)break;
		}
	}
return 0;
}

bessinit(sol,n)
{struct complex z,j,y,h,jprime,yprime,hprime;int i,iv;
z.x=u1;z.y=0.;
/*printf(" bessinit u1=%le u2=%le n=%d\n",u1,u2,n);*/
for(i=0;i<=n;i++)
	{
	bessel(i,&z,&j,&y,&h,&jprime,&yprime,&hprime,&iv);
	if(iv)fprintf(stderr," bessint: pblm %d\n",iv);
	ju1[i]= j.x;
	dju1[i]=jprime.x;
	/*printf(" j,dj[%d]=%le %le\n",i,ju1[i],dju1[i]);*/
	}
z.x=u2;z.y=0.;
for(i=0;i<=n;i++)
	{
	bessel(i,&z,&j,&y,&h,&jprime,&yprime,&hprime,&iv);
	jyu2[i]= (sol==1)? j.x: y.x;
	djyu2[i]= (sol==1)? jprime.x: yprime.x;
	/*printf(" jy,djy[%d]=%le %le\n",i,jyu2[i],djyu2[i]);*/
	}
return 0;
}

int math(xx,qq,r,cv,sol,fnc,norm,f,k) double xx,qq,cv,f[];
int sol,fnc,norm,k[],r;
/* sol= 1 radial,first kind 2 second 3 periodic
at point x, order r, qq, cv=characteristic value a[r](q) or b[r](q)
fnc= 1 b  2 a 3 deriv b 4 deriv a
norm= 1 neutral 2 Ince 3 Stratton
f[3] solution value, largest mag series term last term in summation
k[2] 1 indices of terms in 2nd,3rd element of f */
{
int i,ll,m;
double sqrt(),aux;
if( sol<1 ||sol>3 ||fnc<1 ||fnc>4)goto bomb;
a=cv;qmath=qq; if(qmath!=0.)qinv=1./qmath;
tol=1.e-8;type=  ((fnc%2)<<1) + r%2;
m=coef();
if(m){fprintf(stderr," mathieu:coef failed\n");goto bomb; }
n=r>>1;
p=r%2;
s=mm>>1;
l=ml>>1;
x=xx;
t=1.;
if(sol==3) /* periodic case*/
	{
	switch(fnc)
		{
		case 1: sum(ps);break;
		case 2: sum(pc);break;
		case 3: sum(dps);break;
		default: sum(dpc);
		}
	if(norm==2)
		{
		t= ab[0]; t*=t;
		if(!type)t+=t;
		for(i=0;i<l;i++) {aux=ab[i+1];t+=aux*aux;}
		t=sqrt(t);
		i=m0>>1;
		if(ab[i-1]<0.)t=-t;
		}
	else if(norm>2)
		{
		if(type<=1)
			{
			t=ab[0];
			for (i=0;i<l;i++)
				{
				t+=ab[i+1];
				}
			}
		else
			{
			t=ab[0]*p;
			for(i=0;i<l;i++)
				{
				t+=ab[i+1]*(2*i+p);
				}
			}
		}
	}
else /* radial case*/
	{
	u1= sqrt(qmath)*exp(-x);
	u2=qmath/u1;
	ll=l+s+p;
	bessinit(sol,ll);
	switch (fnc)
		{
		case 1: sum(ds);break;
		case 2: sum(dc);break;
		case 3: sum(dds);break;
		default: sum(ddc);break;
		}
	}
/*printf(" t, total %le %le\n",t,total);*/
f[0]=total/t;f[1]=dmax/t;f[2]=dlast/t;
k[0]=kmax;k[1]=klast;
return m;
bomb:
f[0]=f[1]=f[2]=0.; k[0]=k[1]=0;
return ierrorcode;
}


/* parameter n not used- is this so in original code?*/
int bounds(k,approx,tola,cv,coln) int coln,k;double approx,tola,*cv;
{int m,ka,km;
double dtm,d0,d1,tm;
ka=0;km=k-1;
if(k!=1 && approx<= cv INDEX(0,km-1))a0=1.+cv INDEX(0,km-1);
else
	a0=approx;
while(1)
{/* label 30*/
m=tmofa(a0,&tm,&dtm);
if(m>0)return m;
d0=-tm/dtm;
if(d0==0.)
	{ cv INDEX(0,km)=a0; cv INDEX(1,km)=0.; return -1;}
if(d0>0.)
	{/* a0=lower, search for upper*/
	while(1)
		{
		a1=a0+d0+.1;
		m=tmofa(a1,&tm,&dtm);
		if(m>0)return m;
		d1=-tm/dtm;
		if(d1==0.)
			{cv INDEX(0,km)=a1;cv INDEX(1,km)=0.;return -1;}
		else if(d1>0.)
			{
			a0=a1;d0=d1;ka++;
			if(ka>=4)return 2;
			}   /* inifinite loop if d1<0. a1 doesn't change*/
		else goto l200;
		}
	}
/* a1 upper search for lower*/
while(1)
	{
	a1=a0;d1=d0;
	a0= max(a1+d1-.1,-2.*qmath);
	if(k==1 || a0>cv INDEX(0,km-1))
		{
		m=tmofa(a0,&tm,&dtm);
		if(m>0)return m;
		d0=-tm/dtm;
		if(d0==0.){ cv INDEX(0,km)=a0; cv INDEX(1,km)=0.; return -1;}
		else if(d0>0.)goto l200;
		ka++;
		if(ka>=4)return 2;
		}
	else break;
	}
ka++;
if(ka>=4)return 2;
a0=a1+max(tola,abs(d1));
/* goto 30*/}
l200:a=.5*(a0+d0+a1+d1);
if(a<=a0 || a>=a1)a=.5*(a0+a1);
return m;
}

int mfitr8(tola,cv,dcv) double *dcv,*cv,tola;
{int m,n,last;
double a2,d,tm,dtm;
n=0;last=0;
while(1)
	{
	n++;
	m=tmofa(a,&tm,&dtm);
	if(m>0){*cv=*dcv=0.;return m;}
	d=-tm/dtm;
	if(n>=40 || a-a0 <=tola ||a1-a<=tola ||abs(d)<tola)last=1;
	if(d==0.){*cv=a;*dcv=0.;return m;}
	else if(d<0.)a1=a;
	else a0=a;
	a2=a+d;
	if(last)
		{
		if(a2>a0 && a2<a1)
			{
			m=tmofa(a2,&tm,&dtm);
			if(m>0){*cv=*dcv=0.;return m;}
			*dcv=d=-tm/dtm;
			*cv=a2;
			return m;
			}
		else
			{
			*cv=a;*dcv=d; return m;
			}
		}
	if(a2>a0 && a2<a1) a=a2;
	else a=.5*(a0+a1);
	}
}

/* matheign: compute the eigenvalues
and place in the array eigenv

q= parameter
r= count of desired eigenvalues
odd= 1(true) if odd else 0 for even eigenfunctions

returns the count of eigenvalues actually found
*/
int matheign(q,r,odd,eigenv) int r,odd;double q,eigenv[];
{
int i,j,coln; double *cv;
if(odd)n=r;
else n=r+1;
coln=n;
cv=(double *)malloc( sizeof(double)*coln*6);
j=mfcval(n,r,q,cv,coln);
for(i=0;i<j;i++) eigenv[i]=cv INDEX(0,i) ;
free(cv);
return j;
}


/* return value of mathieu functions
x= argument q=parameter
eigenv= array, previously computed. size r;

sol= 1 radial, first kind
     2 radial, second kind
     3 periodic

fnc= 1 odd solution 
     2 even
     3 derivative of odd
     4 derivative of even

norm=1 neutral
     2 Ince
     3 Stratton
*/

double mathieu(x,q,r,eigenv,sol,fnc,norm)
int sol,fnc,norm,r;double x,q,eigenv[];
{
double fval,f[3];
int k[2],l;
l=math(x,q,r,eigenv[r],sol,fnc,norm,f,k); fval=f[0];
return fval;
}

int mfcval(n,r,qq,cv,coln) double *cv,qq; int n,r,coln;
{int j,k,kk,km,l,m;
double a,dtm,t,tm,tola,sqrt();
tol=1.e-3;
if(n>r) l=2;
else	l=1;
qmath=qq;if(qmath!=0.)qinv=1./qmath;
for(k=1;k<=n;k++)
	{
	j=k;km=k-1;
	if (qmath==0.)
		{
		m=k-l+1;
		cv INDEX(0,km)= cv INDEX(2,km)=cv INDEX(4,km)=m*m;
		cv INDEX(1,km)= cv INDEX(3,km)=cv INDEX(5,km)=0.;
		}
	else if(qmath<0.)
		{ return 0;}
	else
		{
		kk=min(k,4);
		type= ((l%2)<<1)+ (k-l+1)%2 ;
		switch (kk)
			{
			case 1:
				{
				if(qmath<1.)
					{
					if(l==1)
						a=1.-qmath-.125*qmath*qmath;
					else
						{
						a=qmath*qmath;
						a*=(.0546875*a-.5);
						}
					break;
					}
				else
					{
					if(qmath<2.)
						{
						if(l==1)
							a=1.033-1.0746*qmath-.0688*qmath*qmath;
						else
							a=.23-.459*qmath-.191*qmath*qmath;
						break;
						}
					a=-.25-2.*qmath+2.*sqrt(qmath);break;
					}
				}
			case 2:
				{t=l;
				if(qmath*t<6.)
					{
					if(l==1)
						a=4.01521-qmath*(.046+qmath*.0667867);
					else
						a=1.+1.05007*qmath-.180143*qmath*qmath;
					break;
					}
				t=k-1;
				a= cv INDEX(0,km-1)-t+4.*sqrt(qmath);break;
				}
			case 3:
				if(qmath<8.)
					{
					if(l==1)
						a=8.93867+.178156*qmath-qmath*qmath*.0252132;
					else
						a=3.70017+.953485*qmath-.0475065*qmath*qmath;
					break;
					}
				t=k-1;
				a= cv INDEX(0,km-1)-t+4.*sqrt(qmath);break;

			case 4:
				a= cv INDEX(0,km-1)-cv INDEX(0,km-2);
				a=3.*a+cv INDEX(0,km-3);
				break;
			default:
				return (0)/* trouble*/;
			}
		}
	if(qmath<1.)
		{if(k==1)tola= max(min(tol,abs(a)),1.e-7);
		 else tola=tol*abs(a);
		}
	else
		{
		tola=tol*max(qmath,abs(a));
		tola= max(min(.4*sqrt(qmath),min(tola,abs(a))),1.e-7);
		}
	m=bounds(k,a,tola,cv,coln);
	if(m>0)return --j;
	else if(!m)/* m==0*/
		{
		m=mfitr8(tola,&(cv INDEX(0,km)),&(cv INDEX(1,km)));
		if(m>0)return --j;
		}
	t= cv INDEX(0,km)-tola;
	m=tmofa(t,&tm,&dtm);
	if(m>0){cv INDEX(2,km)=cv INDEX(3,km)=0.;continue;}
	else
		{
		cv INDEX(2,km)=t;
		cv INDEX(3,km)=-tm/dtm;
		}
	t= cv INDEX(0,km)+tola;
	m=tmofa(t,&tm,&dtm);
	if(m>0){cv INDEX(4,km)=cv INDEX(5,km)=0.;continue;}
	else
		{
		cv INDEX(4,km)=t;
		cv INDEX(5,km)=-tm/dtm;
		}
	}/* loop over values*/
return j;
}
