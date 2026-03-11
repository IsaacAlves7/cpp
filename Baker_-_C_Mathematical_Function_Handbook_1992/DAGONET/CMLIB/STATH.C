/*
calculator version of stat.c routines
routines for statistical functions and auxilliar functions used in their
computation

chisq() Chi-square
studt() Student's-t distribution
fdist() F distribution
Hot()   Hotelling's T-squared distribution ( generalized t dist)
prob() probability integral
cump() cummulative probability integral (P)
icump() inverse cummulative prob. integral
erf()   error function (related simply to probability integral)
ierf()  inverse of the error function
loggam() logarithm of the absolute value of the gamma function
gamma()  gamma function
beta()    beta function
incbeta() incomplete beta function
incgam() incomplete gamma function
finv()   inverse of F distribution (F value given confidence,d.f.'s)
tinv()   inverse of student's t (t value given confidence,d.f.) two-sided.
chiinv()  "  of chi-sq
ks()	 Kolmogorov-Smirnov distribution (1 sample)
wilcoxon() Wilcoxon cummulative count (integer)
mannw()	 Mann-Whitney confidence (calls wilcoxon)
sign_crit() sign test
invnr()  auxilliary routine to use Newton-Raphson iteration for
		use with finv,tinv.
ncf()   noncentral f distribution
ncc()      "       chi-sq
nct()	   "	   Student's t
incf()	   "        f distribution inverse
inct()     "        t       "         "
incc()     "        chi-sq  "         "

P,Q	PROBABILITY INTEGRALS P+Q=1
L	Bivariate Normal distribution

Copyright 1991 by Louis Baker. All rights reserved.
from More C tools for scientists and engineers by L. Baker
*/
#include <alloc.h>
#include <stdio.h>
#define DOFOR(i,to) for(i=0;i<to;i++)
#define abs(x) ((x)>0.? (x):-(x))
#define max(a,b) ((a)>(b)? (a):(b))
#define min(a,b) ((a)<(b)? (a):(b))
#define pi 3.14159265358979
#define errorcode -1.
#define D sizeof(double)
/*double P(x) double x;{return cump(x);}*/
#define P(x) cump(x)

static double tolrel=1.e-7,tolabsg=1.e-5,tolabsb=1.e-5,
							conabs=1.e-5,conrel=1.e-3;
static int iterkti,iterkt,wilcox;

#define MAXZ 50
main(argc,argv) int argc;char **argv;
{int i,j,k,m,n,df1,df2;
float zz,rr,qq;
double x,Hot(),ierf(),erf(),y,z,q,r,exp(),loggam(),gamma(),ks(),mannw();
double prob(),chisq(),studt(),fdist(),incgam(),incbeta(),cump(),icump();
double finv(),tinv(),chiinv(),ncf(),ncc(),nct(),incf(),inct(),incc(),pow();
double Q(),L(),zn(),Z[ MAXZ ]; int sign_crit();
FILE *fileid,*in;
/*in=fopen("CON:","r");
if(!in)printf(" pblm opening for input\n");*/
/*BEGIN SPECIAL PLOTTING SECTION*/
/*fileid=fopen("PLOT.GAM","w");
fprintf(fileid,"1 \n");
for(i=0;i<100;i++)
	{
	x=(i-50)*.1+.05;
	fprintf(fileid," %f %e\n",x,gamma(x));
	}
fclose(fileid);
fileid=fopen("PLOT.FEW","w");
fprintf(fileid," 6 \n");
for(i=0;i<50;i++)
	{
	x=(i)*.1+.05;
	fprintf(fileid," %f %e %e %e %e %e\n"
,x,loggam(x),erf(x),ierf(x/5.1),prob(x),cump(x),icump(x));
	}
fclose(fileid);
*/
/*
fileid=fopen("PLOT.PRB","w");
fprintf(fileid,"6 \n");
for(i=0;i<50;i++)
	{
	x=(i)*.1+.05;
	y=i*.02;
	fprintf(fileid," %f %e %e %e %e %e %e\n",
x,chisq(x,10.),studt(x,10.),fdist(x,10.,5.),fdist(x,5.,10.)
,incbeta(.5,.5,y),incgam(2.,x));
	}
fclose(fileid);
fileid=fopen("PLOT.inv","w");
fprintf(fileid,"4 \n");
for(i=0;i<50;i++)
	{
	x=(i)*.01+.05;
	fprintf(fileid," %f %e %e %e %e\n"
,x,tinv(x,10.),finv(x,5.,10.),chiinv(x,10.),finv(x,10.,20.));
	}
fclose(fileid);
*/
/*END SPECIAL PLOTTING SECTION.*/
k=1;m=1;q=1.;x=1.;/*defaults*/
iterkt=iterkti=0;
while(1)
{
printf(" enter: 1 for F, 2 chi-sq, 3 Student's 4 Hotelling gen. t\n");
printf(" 5 z->p 6 p->z 7 gamma 8 erf 9 inverse erf 10 prob\n");
printf(" 11 inc beta 12 inc gamma 13 t-inv 14 F-inv 15 chi-sq inv\n");
printf(" 16 noncentral F 17 chi-sq 18 t 19 inv-F 20 inv-t 21 inv-chi\n");
printf(" 22 Kolmogorv-Smirnov 23 Mann Whitney 24 Wilcoxon 25conv param\n");
printf(" 26 P,Q 27 Bivariate Normal L 28 sign test 29 Z\n");
j=scanf("%d", &i);
if(j==-1)printf(" trouble\n");
if(i==2)
	{
	printf(" enter chi-sq value ");j=scanf("%f",&zz);
	z=zz;
	printf(" enter degrees of freedom(integer) ");j=scanf("%d",&k);
	if(k>0)x=k;
	y=chisq(z,x);q=1.-y;
	printf(" z=%f df=%f chisq(z,df)=%e %e\n",z,x,y,q);
	}
else
	if(i==3)
	{
	printf(" enter student t value ");j=scanf("%f",&zz);z=zz;
	printf(" enter degrees of freedom(integer) ");j=scanf("%d",&k);
	if(k>0)x=k;
	y=studt(z,x); q=1.-y;r=q*.5;
	printf(" student's-t(z=%f,df=%f)=%e\n",z,x,y);
	printf("significance levels= %e[1-sided] %e[2-sided]\n" ,r,q);
	r=1.-r;
	printf(" cummulative t[1-sided]=%f %f[two-sided]\n",r,y);
	}
else if(i==1)
	{
	printf(" enter F value ");j=scanf("%f",&zz);x=zz;
	if(j!=-1)printf(" echo F value=%f\n",x);
	else
		printf(" input error\n");
	printf(" enter degrees of freedom(2 integers) ");j=scanf("%d %d",&k,&m);
	if(k>0)z=k;
	if(m>0)q=m;
	y=fdist(x,z,q);
	printf(" F=%f fdist(F,%d,%d)=%e\n",x,k,m,y);
	}
else if(i==4)
	{
		printf(" enter as integers degrees of freedom n and k\n");
		scanf("%d%d",&m,&k);
		if(m>0) n=m;
		printf(" enter value for T (NOT square)\n");
		scanf("%f",&zz);
		x=zz;
		printf(" echo n %d k %d T %f\n",n,k,z);
    	y=Hot(x,n,k);
		printf(" generalized t dist=%e\n",y);
	}
else if (i==5)
	{
	printf(" enter z score\n");
	scanf("%f",&zz);
	x=zz;
	y=cump(x);/*.5*(1.+erf(.7071*x));*/
	printf(" z score=%e\n",y);
	}
else if (i==6)
	{
	printf(" enter p score\n");
	scanf("%f",&zz);
	x=zz;
	y=icump(x);/*1.4142*(ierf(2.*x-1.));*/
	printf(" p score=%e\n",y);
	}
else if (i==7)
	{
	printf(" enter x\n");scanf("%f",&zz);x=zz;
	printf(" gamma(%f)=%f\n", x,gamma(x));
	}
else if (i==8)
	{
	printf(" enter x\n");scanf("%f",&zz);x=zz;
	printf(" erf(%f)=%f\n", x,erf(x) );
	}
else if (i==9)
	{
	printf(" enter x\n");scanf("%f",&zz);x=zz;
	printf(" %f=erf(%f)\n", ierf(x),x );
	}
else if (i==10)
	{
	printf(" enter x\n");scanf("%f",&zz);x=zz;
	printf(" %f=prob(%f)\n", prob(x),x );
	}
else if (i==11)
	{
	printf(" enter x,a,b\n");scanf("%f%f%f",&zz,&rr,&qq);
	x=zz;r=rr;q=qq;
	printf(" %f=inc.beta(%f|%f,%f) %d\n", incbeta(r,q,x),x,r,q,iterkt );
	}
else if (i==12)
	{
	printf(" enter x,a\n");scanf("%f%f",&zz,&qq);
	x=zz;q=qq;
	printf(" %f=inc gamma(%f,%f)\n", incgam(q,x),q,x );
	}
else if (i==13)
	{
	printf(" enter 2-sided confid.,nu\n");scanf("%f%f",&zz,&qq);
	x=zz;q=qq;
	if(x<=0. || x>=1.)
			printf(" confidence must be between 0 and 1\n");
	else
			{
			printf(" %f=inverse t(%f,%f)\n", tinv(x,q),x,q);
			}
	}
else if (i==14)
	{
	printf(" enter confid,d.f.1,2\n");scanf("%f%f%f",&zz,&qq,&rr);
	x=zz;q=qq;r=rr;
	if(x<=0. || x>=1.)
			printf(" confidence must be between 0 and 1\n");
	else
		{
		printf(" %f=inv F(%f,%f,%f)\n", finv(x,q,r),x,q,r);
		}
	}
else if (i==15)
	{
	printf(" enter 1-sided confid.,nu\n");scanf("%f%f",&zz,&qq);
	x=zz;q=qq;
	if(x<=0. || x>=1.)
			printf(" confidence must be between 0 and 1\n");
	else
			{
			printf(" %f=inverse chi-sq(%f,%f)\n", chiinv(x,q),x,q);
			}
	}
else if (i==16)
	{
	printf(" enter  F, df1,df2, noncentrality.df integers\n");
	scanf("%f%d%d%f",&zz,&df1,&df2,&qq);
	x=zz;q=qq;
	if(x<=0. )
			printf(" noncentrality must be positive");
	else
			{
			printf(" noncentral F=%f\n", ncf(x,q,df1,df2) );
			}
	}
else if (i==17)
	{
	printf(" enter  chisq, df,noncentrality.df integer\n");
	scanf("%f%d%f",&zz,&df1,&qq);
	x=zz;q=qq;
	if(x<=0. )
			printf(" noncentrality must be positive");
	else
			{
			printf(" noncentral chisq=%f\n", ncc(x,q,df1) );
			}
	}
else if (i==18)
	{
	printf(" enter  t, df, noncentrality.df integer\n");
	scanf("%f%d%f",&zz,&df1,&qq);
	x=zz;q=qq;
	if(x<=0. )
			printf(" noncentrality must be positive");
	else
			{
			printf(" noncentral t=%f\n", nct(x,q,df1) );
			}
	}
else if (i==19)
	{
	printf(" enter  confidence, df1,df2, noncentrality.df integers\n");
	scanf("%f%d%d%f",&zz,&df1,&df2,&qq);
	x=zz;q=qq;
	if(q<=0. )
			printf(" noncentrality must be positive");
	else
			{
			printf(" inverse noncentral F=%f\n", incf(x,df1,df2,q) );
			}
	}
else if (i==20)
	{
	printf(" enter  confidence, df, noncentrality.df integer\n");
	scanf("%f%d%f",&zz,&df1,&qq);
	x=zz;q=qq;
	if(q<=0. )
			printf(" noncentrality must be positive");
	else
			{
			printf(" inverse noncentral t=%f\n", inct(x,df1,q) );
			}
	}
else if (i==21)
	{
	printf(" enter  confidence, df, noncentrality.df integer\n");
	scanf("%f%d%f",&zz,&df1,&qq);
	x=zz;q=qq;
	if(q<=0. )
			printf(" noncentrality must be positive");
	else
			{
			printf(" inverse noncentral chisq=%f\n", incc(x,df1,q) );
			}
	}
else if(i==22)
	{
	printf(" enter d, n for kolmogorov smirnov\n");
	scanf("%e%d",&zz,&n);z=zz;
	q=ks(z,n);r= (1.-q);
	printf(" answer=%e siglevel=%e\n",q,r);
	}
else if(i==23)
	{
	printf(" enter n,m,u integers for Wilcoxon (Mann-Whitney) 2 sample statistic\n");
	scanf("%d%d%d",&n,&m,&k);
	i= m*n+(m*(m+1)>>1)-k;
	printf(" Wilcoxon=%d T=%d, mann-whitney= %e\n",wilcox,i,mannw(n,m,k));
	}
else if(i==24)
	{
	printf(" enter n,m integers for Wilcoxon matched pair signed rank\n");
	scanf("%d%d",&n,&m);
	j=wilcoxs(n,m);
	printf(" count=%d, prob.= %e\n",j,j/pow(2.,(double)n));
	}
else if(i==25)
	{
	printf(" previous iteration count %d\n",iterkt);
	printf(" enter tolabsg=%e,tolabsb=%e,tolrel=%e nonzero to alter\n"
	,tolabsg,tolabsb,tolrel);
	scanf("%e%e%e",&zz,&qq,&rr);
	if(zz!=0.) {tolabsg=zz;printf(" new tolabsg=%e\n",tolabsg);}
	if(qq!=0.) {tolabsb=qq;printf(" new tolabsb=%e\n",tolabsb);}
	if(rr!=0.) {tolrel=rr;printf(" new tolrel=%e\n",tolrel);}
	printf(" inverse funct. last iteration count=%d\n",iterkti);
	printf("conabs=%e,conrel=%e enter nonzero to alter\n"
	,conabs,conrel);
	scanf(" %e %e",&zz,&qq);
	if(zz!=0.) {conabs=zz;printf(" new conabs=%e\n",conabs);}
	if(qq!=0.) {conrel=qq;printf(" new conrel=%e\n",conrel);}
	}
else if (i==26)
	{
	printf(" enter  x\n");
	scanf("%f",&zz);
	x=zz;
	printf(" P(x)=%le Q(x)=%le\n", P(x),Q(x) );
	}

else if (i==27)
	{
	printf(" enter h,k,rho\n");
	scanf("%f%f%f",&zz,&qq,&rr);
	x=zz;q=qq;r=rr;
	if(r<-1. || r>1.) printf(" correlation must be between -1 and 1\n");
	else
		printf(" L(%le,%le,%le)=%le\n",x,q,r,L(x,q,r));	
	}
else if (i==28)
	{
	printf(" enter n,alpha sign test (two-tailed) critical value\n");
	scanf("%d%f",&n ,&zz);x=zz;
	printf(" %d\n", sign_crit(x,n) );
	}

else if (i==29)
	{
	printf(" enter n,x Z(n)[x]=n-th derivative of Normal prob.\n");
	scanf("%d%f",&n ,&zz);x=zz;if(n>MAXZ)n=MAXZ;
	printf(" %le\n", zn(x,n,Z ));
	/* out of memory=no room to print Z*/
	}


else
	{printf(" terminating j i %d %d\n",j,i);break;}
}
exit(0);
}
double arr(m,n)int m,n;
{/* m+n!/m!n! */
double product;
int i,nn,mm;
if(n>m){nn=n;mm=m;}
else{nn=m;mm=n;}
product=1.;
for(i=1;i<=mm;i++)
	{product*=(((double)(nn+i))/((double)i));
	}
return product;
}
int wilcoxs(ni,mi) int ni,mi;
{
int n,m,x;
if(mi<0)return 0;
if(!mi)return 1;
x=(ni*(ni+1)) >>1;
if(mi>x) return wilcoxs(ni,x);
m=mi;n=ni;
if(n<0 || m<0 )printf(" warn n %d m%d\n",n,m);
return wilcoxs(n-1,m)+wilcoxs(n-1,m-n);
}

int wilcoxon(ni,mi,ui) int ni,mi,ui;
{
int n,m,u;
if(ui>ni*ni) return -1;
if(ui>=mi*ni) return (int)(arr(mi,ni) +.5);
if(ui<0)return 0;
if(!ui)return 1;
if(!mi)return 1;
u=ui;
if(mi>0 && ni>0 && ni<mi) {m=ni;n=mi;}else{m=mi;n=ni;}
if(n<0 || m<0 )printf(" warn n %d m%d\n",n,m);
return wilcoxon(n-1,m,u-m)+wilcoxon(n,m-1,u);
}


double mannw(n,m,u) int n,m,u;
{int a;double arr(),denom;
a=wilcoxon(n,m,u);
wilcox=a;
denom=arr(m,n);
return (double)a/denom;
}

int sign_crit(double two_tail,int n)
{
/* for one-tailed use twice the one-sided value for two_tail*/
double pow(),sum,term,stop;int i;
sum=term=1.;stop =.5*two_tail*pow(2.,(double)n);
for(i=0;i<n;i++)
	{term*= (n-i)/(i+1.);
	sum+=term;
	if(sum>=stop)return i;
	}
return -1;
}




double ks(d,n) double d;int n;
{
int ind,ndt,ndp,ndd,nddp,i,j,k,jmax,na=75;
double dk;
double nd,*q,*f,sum,ft,fu,fv,fn,ci,pow(),sqrt(),exp(),mult;

if(n>na)
	{
	d*=sqrt((double)n);
/*	fv= d+1./(sqrt((double)n)*6.);
	ft= exp(-2.*fv*fv);
	fn=1.-2.*ft;
*/
	fu=exp(-2.*d*d)*(1.-d/sqrt((double)n)*.66666666);
	sum=1.-2.*fu;
/*	printf(" fn, sum %e %e\n",fn,sum);*/
	return sum;
	}
/* na chosen to avoid overflow problems for large n*/
if(n==1)return 2.*d-1.;/* pblms d<.5???*/
nd= n*d;
fn=(double)n;
ndt= 2.*nd;
if(ndt<1)return 0.;
ind=nd;
ndp=ind+1;
ndd= min(n,ind<<1);
f=malloc( D * (n+2) );
q=malloc( D * (n+2) );
nddp=ndd+1;
ci=1.;
f[0]=1.;
for(i=0;i<n;i++)
	{
	f[i+1]=f[i]*ci;
	ci++;
	}
mult=f[n]/pow(fn,fn);
/* might be efficient to invert the f[i] here*/
for(i=0;i<=n;i++) f[i]=1./f[i];
q[0]=1.;
if(ndd)
	{
	ci=1.;
	for(i=1;i<=ndd;i++)
		{
		q[i]=pow(ci,(double)i)*f[i];
		ci++;
		}

	if(ndp>n) goto r1;
	fv= ndp-nd;
	jmax= fv+1;
	for(i=ndp;i<=ndd;i++)
		{
		sum=0.;
		ft=nd;
		k=i;
		fu=fv;
		for(j=0;j<jmax;j++)
			{sum+=pow(ft,(double)j-1)*pow(fu,(double)k)*
				(f[j]*f[k]);
			ft++;fu--;k--;
			}
		q[i]-=2.*nd*sum;
		jmax++;
		fv++;
		}
	if(ndd==n)goto r1;
	}
for(i=nddp;i<=n;i++)
	{
	sum=0.;
	ci=1.;
	ft=2.*nd;
	for(j=1;j<=ndt;j++)
		{
		ft--;
		k=i-j;
		sum+= ci*pow(ft,(double)j)*q[k]*f[j];
		ci=-ci;
		}
	q[i]=sum;
	if(sum<0. || sum*mult>1.)
		{/*trouble*//*printf(" trouble k-s:  sum=%e\n",i,sum);
		for(j=1;j<=i;j++)printf(" q[%d]=%e f=%e\n",j,q[j],f[j]);*/
		return 1.;
		}
	}
r1: 
free(f);free(q);
return q[n]*mult ;
}


double gamma(x) double x;
{
double y,z,exp(),loggam();
y=exp(loggam(x));
if(x>=0.)return(y);
z=  2*(((int)(-x))%2) -1;
return(y*z);
}

double loggam(x) double x;
{
int i;
double z,tmp,ser,log(),sin(),*coeff;
static double logsr2pi=.918938533;
static double b[9]={.035868343,-.193527818,.482199394,-.756704078,
.918206857,-.897056937,.988205891,-.577191652,1.0};

/*if( x<0.&& x> -1. )
	{
	return((loggam(1.+x)-log(-x)));
	}
else requires two levels of recursion and  log call,not sin
*/
if (x<-0.) /*was x< -1. when above implemented */
		{/*transform to x>0. will blow up if x integer, as it should*/
		z=1.-x;/* z>2. */
		return(log(pi/abs(sin(pi*z)))-loggam(z) );
		}
else
	if (x<=1.)/* 0<=x<1 */
		{
		/*z=1.-x*/;/*  0<=z<1*/
		/*return( log(z*pi/sin(pi*z))-loggam(1.+z));*/
		/* Ab& Stegun-takes less than half the time*/
		if(x==0.)return 0.;
		tmp=b[0];
		coeff=&(b[1]);
		for(i=1;i<9;i++)tmp= tmp*x+ *(coeff++);
		return(log(tmp/x));
		}
/* use below for x>1.*/
else
	if(x<=2.)
		{
		tmp=b[0];
		coeff=&(b[1]);
		z=x-1.;
		for(i=1;i<9;i++)tmp= tmp*z+ *(coeff++);
		return(log(tmp));
		}
z=1./x;
tmp=z*z;
/*ser= (1./12.+tmp*(-1./360.+tmp*(1/1260.-tmp/1680.)   ))/x;*/
ser= (.08333333333333+tmp*(tmp*(0.000793650793-.000595238095*tmp)
	-.002777777777))*z;
return (logsr2pi-x+(x-.5)*log(x)+ser);
}
#define small 1.e-30

double incgam(a,x) double x,a;
{
int i,itmax=100;
double gln,exp(),log(),loggam(),sum,ap,del,fi,start,
tol=3.e-7,c0,d0,an,ana,anf,old,offset,mult,delta;
/* error condition return -1 on invalid arguments*/
if( x< 0. || a<0. ) return(-1.);
if(x==0.)return(0.);
gln=loggam(a);
if (x< (a+1.))
	{
	/*series*/
	offset=0.;
	mult=1.;
	ap=a;
	sum=1./a;
	del=sum;
	DOFOR(i,itmax)
		{
		ap++;
		del*=x/ap;
		sum+=del;
		if( abs(del)<abs(sum)*tol) goto fini;
		}
		printf(" trouble incomplete gamma series\n");
	}
else
	{
	offset=1.;
	mult=-1.;
	old=0.;
	start=small;
	sum=start;
	d0=0.;c0=sum;
	DOFOR(i,itmax)
		{
		fi=i;
		if(i)ana=fi;
		else ana=1.;
		d0= (x+d0*ana);
		c0=(x+ana/c0);
		if(d0==0.)d0=small;
		if(c0==0.)c0=small;
		d0=1./d0;
		delta=d0*c0;sum*=delta;
		ana=fi+1.-a;
		d0= (1.+d0*ana);
		c0=(1.+ana/c0);
		if(d0==0.)d0=small;
		if(c0==0.)c0=small;
		d0=1./d0;
		delta=d0*c0;sum*=delta;
		if( abs(delta-1.)<tol)
			{sum-=start;goto fini;}
		}
	printf(" trouble incomplete gamma cont. fract\n");
	}
/*return(-1.);*/
fini:return(offset+mult*sum*exp(-x+a*log(x)-gln));
}

double chisq(csq,nu) double csq,nu;
{
double incgam();
return (incgam(.5*nu,.5*csq));
}


double incbeta(aa,bb,x) double aa,bb,x;
{
int itmax=25,m;
/* uses Abramowitz and Stegun 26.5.8 for Ix(a,b)
26.5.9 seems less reliable */
double offset,bmult,exp(),log(),loggam(),a,b,z;
double error,fm,twicefm,dc,aplusb,am1,ap1, factor,tol=1.e-7,
	h,d,c,delta;
iterkt=0;
if(x<0. || x>1.) return(errorcode);
if(aa==0. )return(1.e10);
if (x==1.||x==0.) return  x ;
bmult= exp(loggam(aa+bb)-loggam(aa)-loggam(bb)+aa*log(x)+bb*log(1.-x));
/*printf(" incbeta x,aa,bb=%e %e %e bt=%e\n",x,aa,bb,bt);*/
 if(x <((aa+1.)/(aa+bb+2.)) )
 	{
 	a=aa;
 	b=bb;
 	z=x;
 	offset=0.;
	bmult/=aa;
	}
else
	{
	a=bb;
	b=aa;
	z=1.-x;
	bmult=-bmult/bb;
	offset=1.;
	};
aplusb=a+b;
am1=a-1.;
ap1=a+1.;
d=0.;h=small;c=h;
for(m=0;m<=itmax;m++)
	{
	fm=(double)(m);
	twicefm=(double)(m<<1);
	if(m)
		dc=fm*(b-fm)*z/((am1+twicefm)*(a+twicefm));/*d2m*/
	else
		dc=1.;
	d=1.+d*dc;
	c=1.+dc/c;
	if(d==0.)d=small;
	if(c==0.)c=small;
	d=1./d;
	delta=d*c;
	h*=delta;
	dc=-(a+fm)*(aplusb+fm)*z/((a+twicefm)*(ap1+twicefm));/*d2m+1*/
	d=1.+d*dc;
	c=1.+dc/c;
	if(d==0.)d=small;
	if(c==0.)c=small;
	d=1./d;
	delta=d*c;
	h*=delta;
	if(abs(delta-1.)<tol)return h*bmult+offset;
	iterkt++;
    };
/*printf(" inc. beta noconv.\n");*/
return (errorcode);/*or return best guess (offset+bmult*conv)*/
}

double beta(a,b) double a,b;
{
double exp(),loggam();
return  exp( loggam(a)+loggam(b)-loggam(a+b));
}

/*26.5.4 series expansion-
use only as a check on incbeta*/
/*
double IncBeta(a,b,x) double a,b,x;
{
int n,nmax=25;
double sum,term,log(),ap1,aplusb,pow,exp(),factor,beta(),fn;
ap1=a+1.;
aplusb=a+b;
factor= exp(a*log(x)+b*log(1.-x))/(beta(a,b)*a);
sum=1.;
pow=x;
iterkt=0;
for(n=1;n<nmax;n++)
	{
	fn=(double)n;
	term= pow*beta(ap1,fn)/beta(aplusb,fn);
	sum+=term;
	if( abs(term/sum)<tolrel || abs(term)<tolabsb)return(sum*factor);
	pow*=x;
	iterkt++;
	}
printf("IncBeta warn\n");
return(sum*factor);
}
*/

double studt(t,nu) double t,nu;
{
double incbeta();
return(1.-incbeta(nu*.5,.5, nu/(nu+t*t)) );
}

double Hot(t,n,k) double t;int n,k;
{
double incbeta(),nu;
nu=n-1;
return(1.-incbeta(nu*.5,.5*k,nu/(nu+t*t)) );
}

double fdist(f,nu1,nu2) double f,nu1,nu2;
{
double incbeta();
return(incbeta(.5*nu2,.5*nu1,nu2/(nu2+nu1*f)) );
}


int g_df1,g_df2;
double g_noncent;

double ncf(f,noncen,df1,df2)
double f,noncen; int df1,df2;
{
/* fdist returns Q want P=1-Q for P' then take 1-P' */
int i,itmax=10;
double nu1,nu2,exp(),fdist(),tol=1.e-3,term,sum,coef,arg1,arg2;
nu1=df1;nu2=df2;
arg1=noncen*.5;
coef=1.;
for (i=0,sum=0.;i<itmax;i++)
	{
	coef /= (double) max(i,1) ;
	term= coef*(1.- fdist(f, nu1+2.*i,nu2) );
	sum=sum+term;
	if( abs(term/sum) < tol)break;
	coef *= arg1;
	}
if(i>itmax)
	{
	printf(" no convergence in noncentral F\n");
	return(-1.);
	}
/*printf(" number of terms for noncentral F=%d\n",i+1);*/
return( 1.- exp(-arg1)*sum);
}

double nct(f,noncen,df1)
double f,noncen; int df1;
{
/* fdist returns Q want P=1-Q for P' then take 1-P' */
int i,itmax=10;
double nu1,nu2,exp(),incbeta(),tol=1.e-3,term,sum,coef,arg1,arg2;
nu1=df1;
arg1=noncen*noncen*.5;
nu2= nu1/(nu1+f*f);
coef=1.;
for (i=0,sum=0.;i<itmax;i++)
	{
	coef /= ((double) max(i*2-1,1) )*((double) max(i*2,1));
	term= coef*(incbeta(nu1*.5,.5+i,nu2) );
	sum=sum+term;
	if( abs(term/sum) < tol)break;
	coef *= arg1;
	}
if(i>itmax)
	{
	printf(" no convergence in noncentral t\n");
	return(-1.);
	}
/*printf(" number of terms for noncentral t=%d\n",i+1);*/
return( 1.- exp(-arg1)*sum);
}

double ncc(f,noncen,df1)
double f,noncen; int df1;
{
/* fdist returns Q want P=1-Q for P' then take 1-P' */
int i,itmax=10;
double nu1,nu2,y,exp(),chisq(),tol=1.e-3,term,sum,coef,arg1,arg2;
nu1=df1;
arg1=noncen*.5;
coef=1.;
for (i=0,sum=0.;i<itmax;i++)
	{
	coef /= (double) max(i,1) ;
	y=chisq(f, nu1+2.*i);
	term= coef*(y );
	sum=sum+term;
	if( abs(term/sum) < tol)break;
	coef *= arg1;
	}
if(i>itmax)
	{
	printf(" no convergence in noncentral chisq\n");
	return(-1.);
	}
/*printf(" number of terms for noncentral chisq=%d\n",i+1);*/
return( exp(-arg1)*sum);
}

double finv(f,nu1,nu2) double f,nu1,nu2;
{
double sqrt(),invnr(),fdist(),icump(),guess,exp(),a,b,h,y,l,w;
	/*guess for F dist*/
if(nu1==1. || nu2==1.)
	{
	if(f>.5)guess=.7;
	else 
	guess= 1./sqrt(f);
	}
else
	{
	y=icump(1.-f);
	l=(y*y-3.)/6.;
	h=2./(1./(nu2-1.)+1./(nu1-1.));
	w= y*sqrt(h+l)/h-(l+5./6.-2./(3.*h))*(1./(nu1-1.)-1./(nu2-1.));	
	guess=exp(w);
/*printf(" debug guess=%f\n");*/
	}
return (invnr(0,f,nu1,nu2,fdist,guess)) ;
}


double ncf3(f,nu1,nu2) double f,nu1,nu2;
{
return ncf(f,g_noncent,g_df1,g_df2);
}

double nct3(f,nu1) double f,nu1;
{
return nct(f,g_noncent,g_df1);
}
double ncc3(f,nu1) double f,nu1;
{
return ncc(f,g_noncent,g_df1);
}

double incf(f,df1,df2,noncen) double f,noncen;
int df1,df2;
{
double nu1,nu2;
double sqrt(),invnr(),ncf(),ncf3(),icump(),guess,exp(),a,b,h,y,l,w;
nu1=df1;nu2=df2;
g_df1=df1;g_df2=df2;g_noncent=noncen;
	/*guess for F dist*/
if(nu1==1. || nu2==1.)
	{
	if(f>.5)guess=.7;
	else 
	guess= 1./sqrt(f);
	}
else
	{
	y=icump(1.-f);
	l=(y*y-3.)/6.;
	h=2./(1./(nu2-1.)+1./(nu1-1.));
	w= y*sqrt(h+l)/h-(l+5./6.-2./(3.*h))*(1./(nu1-1.)-1./(nu2-1.));	
	guess=exp(w);
/*printf(" debug guess=%f\n");*/
	}
return (invnr(0,f,nu1,nu2,ncf3,guess)) ;
}

double tinv(a,nu) double a,nu;
{
double p,t,x;
double y,invnr(),studt(),iguess,icump(),term2,term3,term4,nui;
		/* guess for student's t*/
/*CAVEAT- INAACURATE FOR SMALL NU*/
nui=1./nu;
p=(1.-a)*.5;
t=icump(1.-p);
y=t*t;
term2=t/96.*(3.+y*(16.+5.*y));
term3=t/384.*(-15.+y*(17.+y*(19.+3.*y)));
term4=t/92160.*(-945.+y*(-1920.+y*(1482.+y*(776.+y*79.))));
iguess= t*(1.+(.25*(1.+t*t)+(term2+(term3+term4*nui)*nui)*nui)*nui);
/*attempt itertive improvement*/
return (invnr(1,a,nu,0.0,studt,iguess)) ;
}

double chiinv(a,nu) double a,nu;
{
double p,t,x;
double y,sqrt(),invnr(),chisq(),iguess,icump();
/* approx. for large nu>30 */
x= icump(1.-a);
p=2./(9.*nu);
t= 1.-p+x*sqrt(p);
t=t*t*t;
iguess= nu*t;
return(invnr(1,a,nu,0.0,chisq,iguess));
}

double invnr(narg,x,nu1,nu2,value,iguess)
int narg;
double x,nu1,nu2,(*value)(),iguess;
{
double y,z,zold,znew,dz0=.005,dz,zp,v;
double deriv,delta,deltas,resid,rnew,rold;
double sqrt();
int maxit=40,halvetop=10,i,j;
z=iguess;
iterkti=0;
j=0;
if(narg)
	{ dz0=.05;
	}
for(i=0;i<maxit;i++)
	{
	iterkti++;
	dz= dz0+.001*z;
	zp=z+dz;
/*	v= (narg)? value(z,nu1):value(z,nu1,nu2);
	y= (narg)? value(zp,nu1):value(zp,nu1,nu2);
*/
	if(narg)
		{v=value(z,nu1);
		 y=value(zp,nu1);
		 }
	else
		{v=value(z,nu1,nu2);
		 y=value(zp,nu1,nu2);
		}
/*if(narg)printf(" debug z=%f v=%f y=%f\n",z,v,y);*/
	deriv=(y-v)/dz;
	resid= v-x;	
	delta= -resid/(deriv+.0001);
	resid*=resid;/*square residual as metric*/
	if(resid<conabs|| abs(delta/resid)< conrel) return(z);
/*if(narg)
printf(" debug deriv %f resid %f  delta=%f\n",deriv,resid,delta);*/
	z+=delta;
	deltas=delta;
	z= max(z,0.001);
	z=min(z,1.e4);/*no reasonable F should exceed this*/
	for(j=0;j<halvetop;j++)
		{
/*		v= (narg)? value(z,nu1):value(z,nu1,nu2);*/
		if(narg)
			{v=value(z,nu1);}
		else
			{v=value(z,nu1,nu2);
			}
		rnew=v-x;
		rnew*=rnew;
		if(rnew<resid)break;
/*if(narg)printf(" worse z=%f rn=%f old=%f v=%f %d\n",z,rnew,resid,v,j);*/
			delta*=.5;
			z-=delta;
		if(halvetop-j==1){/* no soap*/
						z-=deltas*.5;/* try other direction*/
						}
		}
	}
		
	printf(" iter=maxit trouble resid=%f z=%f narg=%d,iguess=%f x=%f\n"
	,resid,z,narg,iguess,x);
	if(narg)printf(" nu1=%f\n",nu1);
	else
		printf(" nu1,nu2 %f %f\n",nu1,nu2);
return(z);
}


double erf(x) double x;
{
double exp(),t,z,sign=1.;
if (x<0.)sign=-1.;
x=abs(x);
t= 1./(1.+x*.3275911);
z=((((t*1.061405429-1.453152027)*t+1.421413741)*t-.284496736)*t
	+.254829592)*t;
return (1.-exp(-x*x)*z)*sign;
}

double ierf(y) double y;
{
int i,maxi=20;
double test,sqrt(),pow(),t,z,q,log(),c,x,dx,f,df;
if(y==0.)return(0.);
if(y>=1.)return(1.e10);
if(y<0.) return(-ierf(-y));
c= 1./ sqrt(sqrt(1.-y));
c=sqrt(sqrt(c));/* for higher power version*/
/*guess*/
x=1.;
if(y<1.)x=y*.7;

for(i=0;i<=maxi;i++)
	{
	f=(((((.0000430638*x+.0002765672)*x+.0001520143)*x+
	.0092705272)*x+.0422820123)*x+.0705230784)*x+1.-c;
	df=(((((6.*.0000430638*x+5.*.0002765672)*x+4.*.0001520143)*x+
	3.*.0092705272)*x+2.*.0422820123)*x+.0705230784);


/*lower power version
	f= 1.-c+(((.078108*x+.000972)*x+.230389)*x+.278393)*x;
	df= (((4.*.078108*x+3.*.000972)*x+2.*.230389)*x+.278393);
	*/
	dx= -f/df;
	x=x+dx;
	test=abs(dx);
/*low power
	if( test < 1.e-6 || test< .0001*x)break;
*/
	if( test < 1.e-8 || test< .0000001*x)break;
/*printf(" iter=%d x=%f dx=%f \n",i,x,dx);*/
	}
if(i==maxi)printf(" ierf max it=%d,x=%f,dx=%f,arg=%f\n",i,x,dx,y);
return (x);
}

double prob(x) double x;
{
double erf();
return(erf(x*.707106781));
}

double cump(x) double x;
{/* area under normal curve from -infinty to x*/
static double rt=.7071067812;
double erf();
return(.5*(1.+erf(x*rt)));
}

double icump(x) double x;
{
double ierf();
static double rti=1.414213562;
return(rti*(ierf(2.*x-1.)));
}

double inct(a,df1,noncen) double a,noncen;
int df1;
{
double nu1,nu2;
double invnr(),nct(),icump(),guess;
double nui,t,y,term2,term3,term4,p;
g_noncent=noncen;g_df1=df1;
nu1=df1;nu2=0.;
	/*guess for t dist*/
		/* guess for student's t*/
/*CAVEAT- INAACURATE FOR SMALL NU*/
nui=1./nu1;
p=(1.-a)*.5;
t=icump(1.-p);
y=t*t;
term2=t/96.*(3.+y*(16.+5.*y));
term3=t/384.*(-15.+y*(17.+y*(19.+3.*y)));
term4=t/92160.*(-945.+y*(-1920.+y*(1482.+y*(776.+y*79.))));
guess= t*(1.+(.25*(1.+t*t)+(term2+(term3+term4*nui)*nui)*nui)*nui);
return (invnr(1,a,nu1,nu2,nct3,guess)) ;
}

double incc(a,df1,noncen) double a,noncen;
int df1;
{
double nu1,nu2;
double sqrt(),invnr(),nct(),icump(),guess,exp(),p,t,x;
nu1=df1;nu2=0.;
g_noncent=noncen;g_df1=df1;
/* approx. for large nu>30 */
x= icump(1.-a);
p=2./(9.*nu1);
t= 1.-p+x*sqrt(p);
t=t*t*t;
guess= nu1*t;
return (invnr(1,a,nu1,nu2,ncc3,guess)) ;
}

#define tol 1.e-10

double zn(x,n,z) double x,z[];int n;
{double p;int i;
z[0]=p= exp(-.5*x*x)/sqrt(2.* pi );
if(!n)return p;
z[1]=-p*x;
if(n==1)return z[1];
for(i=2;i<=n;i++)
	{z[i]=-(x*z[i-1]+z[i-2]*(i-1.));
	}
return z[n];
}




double Q(x) double x;
{return 1.-P(x);
}

double L(h,k,rho) double h,k,rho;
{/* L(h,0,rho)*/
int n;
double sum,term,factor,z0[51],z[51];
factor=rho;
term=zn(k,25,z0);
term=zn(h,25,z);
sum=0.;
for(n=0;n<25;n++)
	{
	term=factor*z0[n]*z[n];
	if(abs(term)<tol)break;
	sum+=term;
/*printf(" n=%d z0=%le z=%le term %le sum=%le\n",n,z0[n],z[n],term,sum);*/
	if(abs(term)<abs(sum)*tol)break;
	factor *= rho/(n+1.);
	}
if(abs(term)>tol*abs(sum)){printf(" for L abs. error %le exceeds tol\n",term);}
/*return Q(h)*.5 + sum;L(h,0,rho)*/
return Q(h)*Q(k)+ sum;
}
