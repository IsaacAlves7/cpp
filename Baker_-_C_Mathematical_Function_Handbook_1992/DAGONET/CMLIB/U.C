/*
Confluent Hypergeometric function U for real positive argument
and real parameters

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

Based on N M Temm, Numer. Math. 41,63-82 (1983)

uabx(a,b,x,eps,uprime) returns U and uprime =U'
chu(a,b,x,kmax,eps,u,uprime) where u is an array from 0 to kmax
 giving Gamma(a+k)/Gamma(a) U(a+k,b,x) 0<=k<=kmax and 
 uprime= Gamma(a+kmax)/Gamma(a)U'(a+kmax,b,x)
*/
#include <alloc.h>
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define entier(x) ((x>0.)?(int)x:((int)x)-1)

brec(a,b,k,f,g,x)double a,b,*f,*g,x;int k;
{
int i,km;
double h;km=k-1;
for(i=0;i<km;i++)
	{
	h=*f-*g;
	*g=((i+b)* *g-a* *f)/x;
	*f=h;
	}
return 0;
}

chu(a,b,x,kmax,eps,u,uprime)
double a,b,x,eps,u[],*uprime; int kmax;
{/* gamma(a+k)/gamma(k)*u(a+k,b,x) for k=0,1,kmax and u'(a+k,b,x)*/
int i,j,kk,l,m,p,q,r,s,n,nu,k0,k1,largex; double ar,br,cr,m0,m1;
double pr,qr,rr,sr,vr,bb[9],bb1[9],bx[9],fi[10],sqrt(),sin(),log(),exp();
double c,d,e,f,g,h,w,a1,b1,*v,pow(),delta,t0,t1,gamma(),lou,
u0,u1,u2,u3,x2,y,z,t,mr,er,p2,p1,p0;
if(a<0. || x<0. || kmax<0 || eps<0.){*u=*uprime=errorcode;return;}
if(a==0.)
	{
	u[0]=1;
	*uprime=0.;
	for( i=1;i<=kmax;i++)u[i]=0.;
	return;
	}
if(b<0.)
	{
	v=(double *) malloc( sizeof(double) * (kmax+1));
	c=a-b+1.;
	d=pow(x,-b);
	chu(c,1.-b,x,kmax,eps,v,&w);
	for(j=0;j<=kmax;j++)
		{
		e=(a+j)/d;
		v[j]=-e*v[j];
		d=e/(c+j);
		}
	*uprime=v[kmax];
	u[kmax]=-x*(*uprime*e*w)/(a+kmax);
	c+=x;
	for(j=kmax-1;j>=0;j--)
		{
		u[j]=(-x*v[j+1]+(c+j)*u[j+1])/(a+j);
		}
	return;
	}
if(b>1.)
	{
	n=entier(b);
	b1=b-n;
	a1=a+kmax;
	c=b-a-1.;
	e=c-x;
	m=entier(e);
	if( ((double)m)==e)m--;
	p=    c>=0. && (c==((double)((int)c)));
	q= kmax<=c;
	r=kmax<=m;s=m>=0;
	if(!kmax){r=s=1;}
	if(r){m=kmax; kk= (p?c:n)-m;}
	if(!r)
		{
		if(p&&q)
			{
			g=1;
			i=kmax-1;
			for(j=0;j<=i;j++)
				g*=(j+a);
			f=g*pow(x,-a1);
			g=-a1*f/x;
			brec(a1,a1+1,c-kmax,&f,&g,x);
			}
		else
			{
			chu(a,b1,x,kmax,eps,u,&u3);
			f=u[kmax];
			g=u3;
			brec(a1,b1,n,&f,&g,x);
			if( !p&&s)
				{
				for(j=kmax;j>0;j--)
					u3-=u[j];
					vr=u[0];
					w=u3;
					d=b1;
				}
			n=m+1;
			u[kmax]=f;
			*uprime=g;
			if(!s)n=0;
			}/* else  (if p&&q)*/
		for(j=kmax-1;j>=n;j--)
			{
			h=(-x*g+(j-e)*f)/(a+j);
			g-=f;
			f=u[j]=h;
			}
		}/* !r  */
		if(s)
			{
			if(p)
				{
				vr=pow(x,-a);
				w=-a*vr/x;
				d=a+1.;
				}
			else if(r)
				{
				chu(a,b1,x,0,eps,u,&w);
				vr=u[0];
				d=b1;
				}
			brec(a,d,kk,&vr,&w,x);
			e=b-n-x;
			for(j=0;j<=m;j++)
				{
				if(!j)
					{
					f=u[0]=vr;
					g=w;
					}
				else
					{
					h=-w;
					g=w=-((a+j-1)*vr-(j+e)*w)/x;
					vr=u[j]=h;
					}
				brec(a+j,b+j-m,m-j,&(u[j]),&g,x);
				if(m==kmax)*uprime=g;
				}
			}/* if s */
	}/*b>1*/
else if(x<=1.4 && 0 )
/* defeat with && 0  CAN GIVE poor acc. due to bessel */
	{
	/*printf(" x<1.4 case\n");*/
	n=9;
	vr=12.56637;
	rr=(x-vr*(b+1.))*.5;
	lou=rr*rr+4.*n*x;
	if(lou<0.){fprintf(stderr," U:lou=%e<0\n",lou);lou=0.;}
	d=(vr*n+rr-sqrt(lou))/(2.*n);
	if(d<4.7124){d=4.7124;vr=w=1;}
	else {vr=abs(sin(d));w=pow(vr,-1.-b);}
	w*=exp(.5*x*(1./vr+1./d));
	delta= eps*exp(-.5*x+(n-1.-b)*log(d))/w;
	z=.5/delta;
	vr=.5-b;
	i=n-1;
	for(j=1;j<=i;j++) z*=(j+vr);
	i=0;
	t=sqrt(x)*pow(z,.5/n);
	e=log(delta)+n-n*log(x);
	label: rr=n+b+t;
	sr=1.+b+t;
	pr=log(rr);
	qr=log(sr);
	f=(log(t+.5)-2.*n*log(t)+(rr-.5)*pr-(sr-.5)*qr-e)
		/(1./(t+.5)-2.*n/t+.5*(n-1)/(rr*sr)+pr-qr);
	if(f<0.)
		{
		t-=f;
		t=sqrt(t*t+2.*x);
		i++;
		if(i<10)goto label;
		}
	else
		{k0=1+ entier(t*t/x-a);}
	/*printf("k0=%d t=%e\n",k0,t);*/
	nu= kmax>=k0 ? 1+kmax : k0 ;
	rr=a+nu;
	w=sqrt(x/rr);
	vr=2.*rr*w;
	t0=kbes( vr,-b);
	t1=kbes( vr,1.-b);
	/*printf(" K(vr) for orders -b,1-b,b=%e %e %e %e\n",t0,t1,vr,b);*/
	vr=pow(w,-b);
	bb[0]=bb1[0]=bx[0]=1.;
	u1=fi[0]=vr*t0;
	u0=fi[1]=vr*w*t1;
	x2=x*x;
	bx[1]=-x/12.;
	bx[2]=x2/288;
	bx[3]=-x*(5.*x2-72.)/51840.;
	bx[4]=x2*(5*x2-288.)/2488320.;
	bx[5]=-x*(x2*(7.*x2-1008.)+6912.)/209018880.;
	bx[6]=x2*(x2*(35.*x2-10080.)+279936.)/75246796800.;
	bx[7]=-x*(x2*(x2*(x2*5.-2520.)+176256.)-746496.)/902961561600.;
	bx[8]=x2*(x2*(x2*(x2*5.-4032.)+566784.)-9953280.)/86684309913600.;
	bb[1]=.5;
	bb[2]=(3.*b-1.)/24.;
	bb[3]=b*(b-1.)/48.;
	bb[4]=(b*(b*(b*15.-30.)+5.)+2.)/5760.;
	bb[5]=b*(b*(b*(b*3.-10.)+5.)+2.)/11520.;
	bb[6]=(b*(b*(b*(b*(b*63.-315.)+315.)+91.)-42.)-16.)/2903040.;
	bb[7]=b*(b*(b*(b*(b*(b*9.-63.)+105.)+7.)-42.)-16.)/5806080.;
	bb[8]=(b*(b*(b*(b*(b*(b*(b*135.-1260.)+3150.)-840.)-2345.)-540.)
		+404.)+144.)/1393459200.;
	for(i=1;i<n;i++)
		{
		t0=bb[i];
		t1=bb1[i]=(b-i)*t0;
		for(j=1;j<i;j++)
			{
			t0+=bb[i-j]*bx[j];
			t1+=bb1[i-j]*bx[j];
			}
		t0=bx[i]+b*t0;
		t1+=bx[i];
		fi[i+1]=(x*fi[i-1]+(i-b)*fi[i])/rr;
		u0+=t0*fi[i+1];
		u1+=t1*fi[i];
		/*printf(" i=%d u0 u1 %e %e\n",i,u0,u1);*/
		}
	w=2.*exp(.5*x)/gamma(1.+a);
	u2=w*u0;
	u3=-w*u1;
	vr=a+1.-b+x;
	k1=nu-1;
	for(j=k1;j>0;j--)
		{
		u1=(-x*u3+(vr+j)*u2)/(a+j);
		u3-=u2;
		u2=u1;
		if(j<=kmax) u[j]=a*u2;
		if(j==kmax)*uprime=a*u3;
		}
	u[0]=-x*u3+vr*u2;
	if(!kmax) *uprime=a*(u3-u2);
	/*printf(" answer=%e\n" ,u[0]);*/
	return;
	}

else/* x>1.4 */
	{
	n=a;
	if( (a-(double)n)==0.)
		{
		n--;
		a-=n;
		kmax+=n;
		}
	largex= x>6.5 && a!=b;
	if(largex) mr=1.;
	else
		{
		mr=0.;
		if(a==b)
			{
			m0=a;
			m1=1.;
			}
		else/* gets wrong answers for a!=b, 1.4<x<6.5 must be here!*/
		/* for a approx=b, should bet m0 approx a, m1 approx 1*/
			{
			m0=0.;
			m1=vr=1.;
			for(r=1;vr>eps*m1;r++)
				{
/* Temme vr*=vr/r I believe is typo */
				vr*=x/r;
				m0+=vr;
				vr*=(a+r)/(b+r);
				m1+=vr;
/*printf(" r=%d vr=%e m0 m1 %e %e\n",r,vr,m0,m1);*/
				}
			vr=exp(-x)*gamma(a+1.)/gamma(b+1.);
			m0=vr*(b+a*m0);m1*=vr;
/*			printf(" m0=%e(a=%e)m1=%e(1) %e\n",m0,a,m1,vr);*/
			}
		}/* else*/
	c=a-b;
	cr=2.+c;
	br=x+a+cr;
	p0=0.;
	vr=p1=er=1.;r=0;
	for(ar=a+r;r<=kmax;)
		{
		p2=(br*p1-ar*p0)/cr;
		er*=ar/cr;
		r++;
		if(largex) mr*=(1.+c/r);
		vr=er/p2;
		br+=2.;
		cr+=1.;
		p0=p1;
		p1=p2;
		}
	w=p0*p1/er;

	for(ar=a+r;vr*(w/p0+mr*(2.+a/r))>=eps;ar+=1.)
		{
		p2=(br*p1-ar*p0)/cr;
		er*=ar/cr;
		r++;
		if(largex) mr*=(1.+c/r);
		vr=er/p2;
		br+=2.;
		cr+=1.;
		p0=p1;
		p1=p2;
		}
	c+=1.;
	vr=x+c;
	u2=1.;
	w=0.;
	u3=-2.*r/(x+sqrt(x*(x+4.*r)));
	
	for(r--;r>0;r--)
		{
		if(largex)
			{
			w+=mr*u2;
			mr*=(r+1)/(c+r);
			}

		u1=(-x*u3+(vr+r)*u2)/(a+r);
		u3-=u2;
		u2=u1;
		if(r>=n && r<=kmax) u[r-n]=u2;
		if(r==kmax)*uprime=u3;
		}
	u1=-x*u3+vr*u2;
	u3-=u2;
	vr=a;
	kk=n-1;
	if(kmax==0) *uprime=u3;
	kmax-=n;

	w= largex ? pow(x, -a)/(a*(w+c*u2)+u1):
				pow(x,-b)/(u1*m1-u3*m0);
	for(r=0;r<=kk;r++)vr/=(a+r);
	if(n)kk=0;
	else
		{
		kk=1;
		u[0]=w*u1;
		}
	w*=vr;
	*uprime *=w;
	for(r=kk;r<=kmax;r++) u[r] *=w;
	}
return;
}

double uabx(a,b,x,eps,uprime) double a,b,x,eps,*uprime;
{
double a1,c,p,q,r,u[1]; int j,n;
n= (a<0.)? (int)a -1 :0;
q=a1=a-n;
u[0]=1.;
if(n<0 && a==b)
	{
	if(a1>0.)chu(a1,a1,x,0,eps,u,&q);
	p=u[0];
	r=p-q;
	for(j=1;j<=(-n);j++)
		{
		r*=x;q=(a1-j)*p;
		p=r-q;
		}
	}
else
	{
	/*if(n<0)printf(" a1=%e\n",a1);*/
	if(a1>0.)chu(a1,b,x,0,eps,u,&q);
	/*if(n<0)printf(" q=%e n=%d\n",q,n);*/
	c=1.+a1-b+x;
	a1-=1.;
	p=u[0];
	if(n<0)	for(j=1;j<=(-n);j++)
		{
		/*printf(" j=%d\n",j);*/
		r=(c-j)*p-x*q;
		q=(a1-j)*(q-p);
		p=r;
		}
	}
*uprime=q;
return p;
}

