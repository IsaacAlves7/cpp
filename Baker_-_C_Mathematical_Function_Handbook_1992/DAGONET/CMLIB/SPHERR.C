/*  Spheroidal Wave Functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

#define sign(a,b)   ( (b>0.)? a: -a)
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)<(b)?(b):(a))

static double ratios[100];

extern double leg[100],sphb[100];
/* spherical bessel function i for negative indices n=1 means -1, etc.*/
/* sinb will also work!*/

double simb(x,n)double x;int n;
{/*  use backward in order, forward in -order: -1,-2,-3,etc. function grows*/
int i;double sqrt(),sm1,sm2,mult,sinh(),cosh(),s,c;
s=sinh(x);c=cosh(x);
sphb[0]=s/x;
sphb[1]=sm1=c/x;
if(n==1)return sm1;
sphb[2]=sm2= (s-c/x)/x;
if(n==2)return sm2;
/* forward recursion*/
mult=-3.;
for(i=3;i<=n;i++)
{sphb[i]= sphb[i-2]+sphb[i-1]*mult/x;mult-=2.;}
return sphb[n];
}

/* caveat c is c^2, which is >0 prolate <0 oblate */
/* renamed to avoid name conflict*/
double dssph(c,s) double c,s;
{return s*(s+1.)+c*(2.*s*s*(2.*s+3.)-1.)/((2*s-1.)*(2.*s+1.)*(2.*s+3.));}
double fs(c,s) double c,s;
{return s*(s+1.)+c*(2.*s*s+2.*s-3.)/((2*s-1.)*(2.*s+3.));}
double es(c,s) double c,s;
{double sqrt();
return c*(s+1.)*(s+2.)/((2.*s+3.)*sqrt((2.*s+1.)*(2.*s+5.)));}
double gs(c,s) double c,s;
{double sqrt();
return c/(2.*s+3.)*sqrt(s*(s+1.)*(s+2.)*(s+3.)/((2.*s+1.)*(2.*s+5.)));}
/* based upon Eispack, ste finds eigenvalues of symm. tridi. matrix*/
int ste(d,e,n) double d[],e[];int n;
{int i,l,m,iterkt;double b,c,f,g,p,r,s,eps,sqrt();
eps=1.e-12;
for(l=0;l<n;l++)
	{
	iterkt=0;/*search for small subdiagonal element*/
	iterate:
	for(m=l;m<n;m++)
		{
		if((m==n-1)||abs(e[m])<eps*(abs(d[m])+abs(d[m+1])))break;
		}
	p=d[l];
	if(m!=l)
		{
		if(iterkt>=30){return -n-l;/*error*/}
		iterkt++;
		/*shift*/
		g=(d[l+1]-p)*.5/e[l];
		r=sqrt(g*g+1.);
		g=d[m]-p+e[l]/(g+sign(r,g));
		s=c=1.;p=0.;
		for(i=m-1;i>=l;i--)
			{
			f=s*e[i];b=c*e[i];
			if(abs(f)>=abs(g))
				{
				c=g/f;
				r=sqrt(c*c+1.);
				e[i+1]=f*r;
				s=1./r;
				c*=s;
				}
			else
				{
				s=f/g;
				r=sqrt(s*s+1.);
				e[i+1]=g*r;
				c=1./r;
				s*=c;
				}
			g=d[i+1]-p;
			r=(d[i]-g)*s+2.*c*b;
			p=s*r;
			d[i+1]=g+p;
			g=c*r-b;
			}
		d[l]-=p;
		e[l]=g;
		e[m]=0.;
		goto iterate;
		}
	else
		{
		if(!l){d[0]=p;}
		else
			{for(i=l;i>=1;i--)
				{
				if(p>=d[i-1]){d[i]=p;break;}
				d[i]=d[i-1];
				}
			}
		}
	}
return 0;
}

double diag[100],sdiag[100];
double ab[100],bl[100];
/* solve, figi set up eigenvalue problem for spheroidal wave functions*/
int solve(c,prolate,odd,order, nmax)double c;int prolate,odd,nmax,order;
{
double cc, (*dia)(),(*subd)();
/* odd=0 even 1 odd*/
int i,offset,j;double s;
if(nmax>100)
	{fprintf(stderr," nmax too big\n");return -1;}
if(((double)nmax)*((double)nmax)< 2.*c)
	{fprintf(stderr," nmax too small for accuracy\n");}
offset=0;
if(!prolate)cc=-c;
else	cc=c;
if(order<2)
	{
	if(order && !odd)offset=2;
	if(order){dia=fs;subd=gs;}
	else {dia=dssph;subd=es;}
	for(i=0;i<nmax;i++)
		{
		j=2*i+odd+offset;s=j;
		diag[i]= dia(cc,s);
		sdiag[i]=subd(cc,s);
		}
	}
else
	{
	for(i=0;i<nmax;i++)
		{
		j=2*i+odd;
		bl[i]= cc*j*(j-1)/((2*(order+j)-3)*(2*(order+j)-1));
		ab[i]= cc*(2*order+j+2)*(2*order+j+1)/((2*(order+j)+3)*(2*(order+j)+5));
		diag[i]=(j+order)*(j+order+1)
		+cc*(2*(order+j)*(order+j+1)-2*order*order-1.)
		/((2.*(order+j)-1.)*(2.*(order+j)+3));
		}
	i=figi(ab,bl,nmax);
	if(i)fprintf(stderr," figi pblm %d\n",i);
	}
i=ste(diag,sdiag,nmax);
if(i)fprintf(stderr," ste pblme index %d",-(i+nmax));
return i;
}

int figi(above,below,nmax)int nmax;double above[],below[];
{
int i,j;   double p,sqrt();
/* note that first nonzero below diagonal element is below[1],not below[0]
 and that there are nmax-1 super/subdiagonal elements each, not nmax*/
for(i=0;i<nmax-1;i++)
	{j=i+1;
	p=above[i]*below[j];
	if(p<0.)
	   {fprintf(stderr," swf:figi: i %d %e %e\n",i,above[i],below[j]);
		return -nmax-i;}
	else if(p==0.)
		{if( above[i]!=0. || below[j]!=0.)
			{
		fprintf(stderr," swf|figi: i %d %e %e\n",i,above[i],below[j]);
			return -3*nmax-i;}
		}
	sdiag[i]=sqrt(p);
	}
return 0;
}

static int nmmodd,firstrho;
double lmn,csquare;
double nsph[100],dsph[100],dsphn[200];

double nmn(m,n,limit) int m,n,limit;
{int odd,k;double r,sum,factor,gamma(),d,twom;
odd=(n-m)%2;r=odd;/*r= 0 or 1.*/
twom=(m<<1);
factor=gamma(twom+r+1.)/gamma(twom+2.*(r+1.));
for(sum=0.,k=0;k<limit;k++,r+=2.)
	{d=dsph[k];
	sum+=d*d*factor;
	factor*=(twom+r+1.)*(twom+r+2.)/((twom+2.*r+3.)*(twom+2.*r+5.));
/*printf(" nmm sum factor now %e %e\n",sum,factor);*/
	}
return 2.*sum;
}

double rho(m,n,limit) int m,n,limit;
{int odd,k;double c,r,sum,factor,gamma(),d,twom,pow(),sqrt();
odd=(n-m)%2;r=odd;/*r= 0 or 1.*/
twom=(m<<1);
c=sqrt(abs(csquare));
factor=gamma( twom+r+1.)/gamma((r+1.));
for(sum=0.,k=0;k<limit;k++,r+=2.)
	{d=dsph[k];
	sum+=d*factor;
	factor*=(twom+r+1.)*(twom+r+2.)/((r+1.)*(r+2.));
/*printf(" rho sum factor now %e %e\n",sum,factor);*/
	}
return ((m-n)%2?-1.:1.)* pow(c,m)/sum;
}


double gamsph(m,r) int m,r;
{
double f,g;
f=m+r;g=2.*f;
return f*(f+1.)+.5*csquare*(1.-(4*m*m-1.)/((g-1.)*(g+3.)));
}

double betsph(m,r) int m,r;
{
double f,g,h;
g=2*(m+r); h=2*m+r;
f=g-1.;f*=f;
return r*(r-1)*(h)*(h-1.)*csquare*csquare/((g-3.)*(g+1.)*f);
}

double asph(m,r) int m,r;
{double twom,g,h;
twom=2.*m;g=twom+r;h=twom+2.*r;
return csquare*g*(g-1.)/((h-1.)*(h+1.));
}
double csph(m,r) int m,r;
{double twom,h;
twom=2.*m;h=twom+2.*r;
return csquare*(r+1.)*(r+2.)/((h+3.)*(h+1.));
}
double bsph(m,r) int m,r;
{double twom,h,f;
twom=2.*m;h=twom+2.*r;f=m+r;
return f*(f+1.)-lmn+csquare*(2.*f*(f+1.)-twom*m-1.)/((h-1.)*(h+3.));
}

int tridi(bl,diag,ab,c,n) int n; double bl[],diag[],ab[],c[];
{  int i;double pivot;
if(n==1)
	{c[0]/=diag[0];return 0;}
for(i=1;i<n;i++)
	{
	pivot=diag[i-1];
	if(pivot==0.){fprintf(stderr," tridi:bad pivot\n");return -1;}
	pivot= bl[i]/pivot;
	diag[i]-=pivot*ab[i-1];c[i]-=pivot*c[i-1];
	}
c[n-1]/=diag[n-1];
for(i=n-2;i>=0;i--)
	c[i]= (c[i]-ab[i]*c[i+1])/diag[i];
return 0;
}

/*double nfwd(m,r) int r,m;
{
if(r<2){return -1;}
if(r==2) return gamsph(m,0)-lmn;
if(r==3) return gamsph(m,1)-lmn;
return gamsph(m,r)-lmn-betsph(m,r)/nfwd(m,r-2);
}                 */

#define tol 1.e-10

int setd(mm,n)int mm,n;
{/* starting from n<(l/2)<100 set the d coefficients by backward recursion*/
int i,j,k,r,s,t,l,m; double sum,alpha,gamma(),temp,scale,pow(),legendrea();
double p,pold,polder,q,qold,qolder,a,b,ratio,newrat;
l=32;
k=l;
m=abs(mm);/* d's are same for pos, neg m except for a scaling done in angular()*/
if(n<m){fprintf(stderr," swf:setd: n=%d<m=%d\n",n,m); return -1;}
nmmodd=(n-m)%2 ;/*? 1: 0; don't need the ?1:0  test*/
nsph[0]=1.;
if(k+1>=100){fprintf(stderr," swf:setd:l too big in setd\n");return -1;}
/* fill d array with N */
nsph[k+1]=dsph[k+1]=0.;
/* down to i=2 r=2i+nmmodd=4 or 5 */
for(i=k;i>=1;i--)
	{r=(i<<1)+nmmodd;nsph[i]= betsph(m,r)/(gamsph(m,r)-lmn-nsph[i+1]);
	}
sum=nsph[1];
nsph[1]=gamsph(m, nmmodd)-lmn;/* special treatment not needed? */
if( abs(sum-nsph[1])> 1.e-4)
	{fprintf(stderr," swf:setd:trouble nsph[1]=%e %e\n",sum,nsph[1]);
	return -1;}
/* temporary norm of d s.t. d[0]=1*/
dsph[0]=1.;
j=(m<<1);
for(i=1;i<=k;i++)
	{
	r=nmmodd+(i<<1); s=j+r;t=j+(r<<1);
	dsph[i]=-dsph[i-1]*nsph[i]*((double)(t-1)*(t+1))/(csquare*(s)*(s-1));
	if(abs(dsph[i])<tol){k=i;break;}
	}
sum=0.;
/* flammer normalization corrected sign of N*/
if(!nmmodd)
	{
	alpha= gamma(2.*m+1.)/gamma((double)m+1.);
	for (i=0;i<k;i++)
		{
		r=i<<1;
		temp=alpha;
		sum+=temp;
		if(abs(temp)<tol*abs(sum)){k=i;break;}
		scale=-.25*(double)((r+2+j)*(r+1+j))/((double)((i+1)*(i+1+m)));
		alpha=temp*scale*dsph[i+1]/dsph[i];
		}
	alpha= gamma((double)n+m+1)*(((n-m)>>1)%2?-1:1)
	/(sum*pow(2.,(double)n-m)*gamma(.5*(n-m)+1.)*gamma(.5*(n+m)+1.));
	}
else
	{
	alpha=.5*gamma(2.*m+3.)/gamma(m+2.);
	for (i=0;i<k;i++)
		{
		r=1+(i<<1);
		temp=alpha;
		sum+=temp;
		if(abs(temp)<tol*abs(sum)){k=i;break;}
		scale=-.25*(r+2+j)*(r+3+j)/((i+1)*(i+2+m));
		alpha=temp*scale*dsph[i+1]/dsph[i];
		}
	alpha= gamma((double)n+m+2)*(((n-m-1)>>1)%2?-1:1)
	/(sum*pow(2.,(double)n-m)*gamma(.5*(n-m-1)+1.)*gamma(.5*(n+m+1)+1.));
	}
for(i=0;i<=k;i++) dsph[i]*=(alpha);
/* set drm for r=-2m or -2m+1*/
if(!m){dsphn[0]=dsph[0];}
else
	{/* setup tridiag  bl  diag ab   = sdiag*/
	i=0;bl[0]=0.;
	r=nmmodd-2;/*r= -2 or -1*/
	diag[0]=bsph(m,r); ab[0]=csph(m,r-2);dsphn[0]=-dsph[0]*asph(m,r+2);
/*printf(" r= i=%d %d %e %e %e %e\n",r,i,dsphn[i],diag[i],ab[i],bl[i]);*/
	i=1;/* r initially -4 or -3 in for loop*/
	for(r=nmmodd-4;r>=nmmodd-j;r-=2,i++)
		{
		dsphn[i]=0.;
		diag[i]=bsph(m,r);
		ab[i]=csph(m,r-2);
		bl[i]=asph(m,r+2);
/*printf(" r= i=%d %d %e %e %e %e\n",r,i,dsphn[i],diag[i],ab[i],bl[i]);*/
		}
/*printf(" calling tridi with i=%d\n",i);*/
	tridi(bl,diag,ab,dsphn,i);
/*for(l=0;l<i;l++)printf(" dsphn[%d]=%e\n",l,dsphn[l]);*/
	}
/* set d rho|r. first value */
if(nmmodd)
	sum=1./((1.-j)*(3.-j));
else
	sum=1./((1.-j)*(1.+j));
pold=0.;polder=1.;qolder=0.;qold=1.;ratio=0.;
for(i=m+1;i<100;i++)
	{
	r=nmmodd-(i<<1);
	if(i==(m+1))a=1.;
	else a= -csph(m,r)*asph(m,r+2);
	b=bsph(m,r);
	p=a*polder+b*pold;
	q=a*qolder+b*qold;
	newrat=p/q;
	if(abs(ratio-newrat)< 1.e-6)break;
	ratio=newrat;
	polder=pold;pold=p;qolder=qold;qold=q;
	}
r=m;if(!m)r=1;
dsphn[r]=-dsphn[r-1]*csquare*sum*newrat;
s= m? 0:1;
firstrho=r;
/*printf(" FIRST d-rho: dsphn[%d]=%e \n",r,dsphn[r]);*/
for(i=m+2;i<100;i++)
	{
        pold=0.;polder=1.;qolder=0.;qold=1.;ratio=0.;
	for(l=i-1;l<25;l++)
		{
		r=-((l<<1)+2-nmmodd);
		if(l==i-1)a=asph(m,r+2);
		else a= -csph(m,r)*asph(m,r+2);
		b=bsph(m,r);
		p=a*polder+b*pold;
		q=a*qolder+b*qold;
		newrat=p/q;
		if(abs(ratio-newrat)< 1.e-6)break;
		ratio=newrat;
		polder=pold;pold=p;qolder=qold;qold=q;
		}
	dsphn[i-1+s]=-dsphn[i-2+s]*newrat;
	if(abs(dsphn[i-1])<1.e-10)break;
	}
return k;
}

static int nlimit;

double angular(mm,n,kind,eta,imag)double eta;int kind,mm,n,imag;
{int i,r,btm,nn,m;
double sum,qlm(),qli(),qbig,qbigger,q;
/*first set up d array*/
/*if(abs(eta)<=1.)*/
{
m=abs(mm);
nlimit=setd(m,n);
btm= (n-m)%2;
sum=0.;
if(kind==1)
	{
	if(imag){for(i=0;i<(nlimit<<1)+m;i++)
			leg[i]=legendrea(m,i+m,eta,0);
		}
	else legptable(eta,(nlimit<<1)+m,m);
	for(i=0,r=btm;i<nlimit;i++)
		{sum+=dsph[i]*leg[r];r+=2;}
	if(m!=mm)
		{/*negative m*/
		sum*=(m%2?-1.:1.)*gamma(n-m+1.)/gamma(n+m+1.);
		}
	return sum;
	}
/* kind==2*/
if(abs(eta)==1. && !imag)return errorcode;
if(!imag)leg2(eta,m,nlimit,6,leg);
else {/*printf("Oblate Radial case: imag. eta for leg2 %e\n",eta);*/
	/* need Q(ix) here */
	if(m || eta!=0.)qleg(m,(nlimit<<1)+m,eta,0,ratios,leg);
	else for(i=0;i<=(nlimit<<1)+m;i++)leg[i]=qli(i,eta);
     }
for(i=0,r=m+btm;i<nlimit;i++)/* r>=0 Q terms*/
	{
	sum+=dsph[i]*leg[r];r+=2;
/*printf(" q sum now %e latest factors %e %e\n",sum,dsph[i],leg[r-2]);*/
	}
if(m)/* skip if m==0*/
	{
	if(m==1)
		{
		if(btm){q=leg[0];qbig=leg[1];qbigger=leg[2];}
		else/*m=1, (n-m)even, dsphn[0]*Q-1,m */
			{
			qbig=leg[0];qbigger=leg[1];
			q=(eta*qbig-(1-m)*qbigger)/m;/* n=-1*/
			}
		}
	else/* for m>1, dsphn[0] multiples Qn,m n>=0 */
		{
		r=btm-(m%2);/*r=-1 or 0*/
		if(!r)r=-2;/* r = -1 or -2 highest neg. index*/
		nn=m+r;/* nn for Qnn,m of highest neg. dsphn[0] term*/
		q=leg[nn];qbig=leg[nn+1];qbigger=0.;
/*printf(" q,qbig,qbigger %e %e %e %d %d nn=%d\n",q,qbig,qbigger,r,btm,nn);*/
		}
	for(i=0;i<firstrho;i++)
		{sum+=dsphn[i]*q;
		nn--;
		if(!(m+1+nn))break;
		qbigger=qbig;
		qbig=q;
		q=(eta*(2*nn+3)*qbig-qbigger*(2+nn-m))/(m+nn+1);
		nn--;
		if(!(m+1+nn))break;
		qbigger=qbig;
		qbig=q;
		q=(eta*(2*nn+3)*qbig-qbigger*(2+nn-m))/(m+nn+1);
		}
	}
/*if(imag)printf(" before p terms %d %d sum=%e\n",nlimit,m,sum);*/
/* need P(ix) for oblate radial case*/
if(imag){for(i=0;i<(nlimit<<1)+m;i++) leg[i]=legendrea(m,i+m,eta,0);}
else legptable(eta,(nlimit<<1)+m,m);
r=(1-btm);
/*if(imag)printf(" legp/legendrea returned r=%d firstrho=%d\n",r,firstrho);*/
for(i=firstrho;i<nlimit;i++)
	{ sum+=q=dsphn[i]*leg[r];r+=2;
/*printf(" P sum now %e term=%e P=%e dsphn=%e i=%d r=%d\n",sum,q,leg[r-2],dsphn[i],i,r-2);*/
	}
return sum;
}
}

double radswf(m,n,kind,eta)double eta;int kind,m,n;
{int i,r,btm,twom,k,limit;
double sum,qlm(),etasq,coef,c,factor,fsave,pow(),gamma(),z;
/*first set up d array*/
if( (abs(eta)>1.&&csquare>0.) ||(csquare<0.))
{
/*printf(" eta %e csquare %e m n kind %d %d %d radial\n",eta,csquare,m,n,kind);*/
twom=m<<1;
nlimit=setd(m,n);
btm= (n-m)%2;
sum=0.;
if(!m) factor=1.;
else
	{if(btm)factor=gamma(twom+2.);
	 else factor=gamma(twom+1.);}
r=btm;fsave=factor;
for(k=0;k<(nlimit);k++)
	{sum+= dsph[k]*factor;
	r+=2.;       factor*=(twom+r)*(twom+r-1.)/(r*(r-1.));
	}
etasq=eta*eta;
if(csquare>=0.){ coef=1.;c=sqrt(csquare);}/*prolate*/
else {coef=-1.;c=sqrt(-csquare);}/*oblate*/
coef=(etasq-coef)/etasq;
coef= pow(coef,.5*m)/sum;
sum=0.;
limit=(nlimit<<1)+m;
z=c*eta;
if(kind==1)
	{
	if(csquare>0.)sjn(z,limit);
	else {sinb(z,limit);for(k=1;k<=limit;k+=2)sphb[k]=-sphb[k];}
	}
else
	{
	printf(" CAVEAT:poor convergence 2nd kind RADIAL SPHER.WAVE FUNCT\n");
	if(csquare>0.)syn(z,limit);
	else {sinb(z,limit);return errorcode;}
	/* oblate radial,second kind not supported through radial()*/
	}
factor=fsave;
for(i=0,r=btm;i<nlimit;i++)
	{sum+=factor*dsph[i]*sphb[r+m];r+=2;
/*printf(" r=%d m=%d %d bes=%e sum %e \n",r,m,r+m,sphb[r+m],sum);*/
	  factor*=-(twom+r)*(twom+r-1.)/(r*(r-1.));}
return sum*coef;
}
/* eta<1.*/
return errorcode;
}

double radial(m,n,kind,xi) int m,n,kind; double xi;
{double coef,angular(),sphjoin(),z;    int limit,imag;
limit=setd(m,n);
coef=sphjoin(m,n,kind,limit);
if(csquare<0.)imag=1;
else imag=0.;
if(coef==0.){fprintf(stderr," warn sphjoin=0. in radial\n");
	return errorcode;}
z=angular(m,n,kind,xi,imag)/coef;
/*printf("radial =%e kappa=%e\n",z,coef);*/
return z;
}

double sphjoin(m,n,kind,limit) int m,n,kind,limit;
{
double sign,sqrt(),c,coef,pow(),r,gamma(),sum,factor,twom;int k,btm,oblate,ms;
/*printf(" in sphjoin %d %d %d %d\n",m,n,kind,limit);*/
/* including sign of factors of -i^ -m,-m+1 (1 even/odd) 1-m,2-m oblate
but NOT factor of i if imaginary*/
btm=(n-m)%2;twom=m<<1;
sum=0.;
if(!m) factor=1.;
else
	{if(btm)factor=gamma(twom+2.);
	 else factor=gamma(twom+1.);}
r=btm;
for(k=0;k<(limit);k++)
	{sum+= dsph[k]*factor;
	r+=2.;       factor*=(twom+r)*(twom+r-1.)/(r*(r-1.));
	}
if(csquare<0.){c=sqrt(-csquare);oblate=1;ms=m%4;}
else {c=sqrt(csquare);oblate=0;}
if(kind==1)
	{
	coef=1.;
	if(btm)
		{
		if(oblate) coef= (ms==1 || ms==2)?-1.:1.;
		return coef*sum*gamma(2.+n+m)*(twom+3.)
		/(pow(2.,(double)(n+m))*dsph[0]*pow(c,1.+m)*gamma(1.+m)
		*gamma(1.+.5*(n+m+1))*gamma(1.+.5*(n-m-1)));}
	{
	if(oblate) coef= (ms)>1?-1.:1.;
	return coef*sum*gamma(1.+n+m)*(twom+1.)
	/(pow(2.,(double)(n+m))*dsph[0]*pow(c,(double)m)*gamma(1.+m)
	*gamma(.5*(n+m)+1.)*gamma(.5*(n-m)+1.));}
	}
/* kind==2*/
	coef= firstrho-1 >=0 ? dsphn[firstrho-1] : dsphn[0];
	sign=1.;
/*printf(" coef=%e sum=%e\n",coef,sum);*/
	if(btm)/*odd*/
		{if(oblate)sign=((ms)<2?-1.:1.);
		factor= -sum*pow(2.,(double)n-m)*gamma(twom+1.)
		*gamma(.5*(n-m-1)+1.)*gamma(.5*(n+m+1)+1.)*coef
		/((twom-3.)*(twom-1.)*gamma(m+1.)*gamma(2.+n+m)*pow(c,m-2.));
		return factor*sign;}
		{if(oblate)sign=(ms==3||!ms?-1.:1.);
		return sum*pow(2.,(double)n-m)*gamma(twom+1.)
		*gamma(.5*(n-m)+1.)*gamma(.5*(n+m)+1.)*coef*sign
		/((twom-1.)*gamma(m+1.)*gamma(1.+n+m)*pow(c,m-1.));}

}
