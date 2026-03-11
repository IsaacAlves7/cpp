
/* functions for the evaluation of the
	exponential integral and the sine and cosine integrals and 
	their relatives

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

e1		E1(x) rational approx.
		error:  < 2x10^-7 0<x<1, < 2x10^-8 1<x<infinity
si,ci	sine and cosine integrals, rational approx.
		error:  < 5x10^-7 0<x<inf
f,g	auxilliary functions for above: not intened for direct use. static.
en		En exponential integral-continued fraction used. aux. routine
e1s	E1, series
eicf	Ei(x), continued fraction approx. aux. routine
e		E(x,n) Exponential integral-general
ei		Ei(x)  Ei-general
Eict 	Ei, moderate and large x, Cody & Thatcher rational approx.
Eis	Ei, small x, Cody & Thatcher rational approx.
eias	asymptotic approximate fit aux. routine
li		Li(x)
alpha	alpha(x,n)
beta	beta(x,n)  Do not confuse with Beta function!

iterated sine and cosine integrals from Ch12 of Bowman et al.:
f(),g() are those of Abramowitz & Stegun Ch. 5
c() and s() of Bowman et al. are g, -f respectively.
Note the following relations:
Abramowitz and Stegun:
Si= integral[0 to x] sin(t)/t
Ci= gamma+log(x)+ integral[0 to x] (cos(t)-1)/t
si=Si- pi/2
f= Ci sin - si cos
g= -Ci cos - si sin

Bowman et. al./Einarsson: rational approx. iterated Si,Ci integrals.
T= c-is   (= g+if in A&S notation)
 = integral[x to infinity] exp i(t-x) /t = -exp ix [Ci+Si-i pi/2]

l2(x,y)	= integral[0 to infinity] T(y+ty) exp ixt/(1+t)
	= c2-is2
	errors: < 2.e-5 for 2<x, 0<=y<1
with l2(x,x)= c01-is01
	errors: <  8.e-10 for 0<x<=2, 9.e-9 x>2

cei	complex exponential integral (all double arg)
ceii	scaled complex exponential integral (all double arg)	
cexpint complex exponential integral (complex arg)
*/

#include "cmlib.h"
#include "protom.h"
#include <stdio.h>
#define DOFOR(i,to) for(i=0;i<to;i++)

/* make external global double for use elsewhere*/


double e1(x) double x;
{
double log(),exp(),y;
	if(x< 0.)return errorcode;
if(x<1.)
	y=(((((.00107857*x-.00976004)*x+.0551968)*x-.24991055)*x
       +.99999193)*x-Egamma)-log(x);
else
y=(((((x+8.5733287401)*x+18.0590169730)*x+8.6347608925)*x
   +.2677737343))/(((((x+9.5733223454)*x+25.6329561486)*x
     +21.0996530827)*x+3.9584969228)*x)*exp(-x);
     return(y);
}

double si(x) double x;
{
double z,cos(),sin(),f(),g();
if(x<1.)
	{
	z=x*x;
/*	return(x*(1.-z*(.055555555-z*(.00166666666-z*.0000283446))) );*/
	return(x*(.999999998-z*(.055555480-z*(.001666289-z*.000027739))) );
	}
else
	return(pi*.5-f(x)*cos(x)-g(x)*sin(x));
/*next for the DeSmet C compiler-others will warn of unreachable stmt*/
return(0.);
}

/*Warning: Name conflict with console input function ci() in DeSmet C*/
double ci(x) double x;
{
double z,y,f(),g(),log(),cos(),sin();
if(x<1.)
	{
	z=x*x;
/*	y=-z*(.25+z*(z*(.0002314814-z*.0000031001)-.01041666))*/
	y=-z*(.25+z*(z*(.000231447-z*.000003046)-.01041666))
     + Egamma +log(x);
	}
else
	y=f(x)*sin(x)-g(x)*cos(x);
return(y);	
}	
/* f, g not for direct user use */

static double f(x) double x;
{ double z;
z=x*x;
return(((((z+38.027264)*z+265.187033)*z+335.67720)*z
     +38.102495))/(((((z+40.021433)*z+322.624911)*z+
     570.236280)*z+157.105423)*x);
}

static double  g(x) double x;
{
double z;
z=x*x;
	return(((((z+42.242855)*z+302.757865)*z+352.018498)*z+
       21.821899))/(((((z+48.196927)*z+482.485984)*z+1114.978885
      )*z+449.690326)*z);
}

#define small 1.e-20

int itktg; /* global to report iteration count*/

double en(x,n) double x,n;
{
int i,itmax=50;
double d,c,delta,h,a,tol=1.e-7,exp();
/* continued fraction: not reliable for small x*/
if(x==0.)
	{
	if(n==1.)return(-1.);/*error*/
	return 1./(n-1.);
	}
h=small;d=0.;c=h;
for (i=0;i<itmax;i++)
	{
	if(!i)a=1.;
	else
		a=i;
	d=x+a*d;
	c=x+a/c;
	if(d==0.)d=small;
	if(c==0.)c=small;
	d=1./d;
	delta=d*c;
	h*=delta;
	a= n + i;
	d=1.+a*d;
	c=1.+a/c;
	if(d==0.)d=small;
	if(c==0.)c=small;
	d=1./d;
	delta=d*c;
	h*=delta;	
	if(abs(delta-1.)<tol)break;
	}
return(h*exp(-x));/* return best guess*/
}

double e1s(x) double x;
{/* e1, series*/
int n; double sum,factor,m,tol=1.e-8;
sum=0.;factor=-x;
for(n=1;n<1000;n++)
	{
	sum+=factor;
	if( abs(factor)<tol*abs(sum))break;
	m=n+1.;
	factor *= -x*n/(m*m);
	}
return  -Egamma-log(x)-sum;
}

double eicf(x) double x;
{
int i,itmax=100;
double podd,peven,qodd,qeven,y,z,a,b,tol=1.e-6,zold,exp();
qodd=0.;qeven=1.;
podd=1.;peven=0.;
z=0.;
/* continued fraction not reliable for small x*/
for (i=0;i<itmax;i++)
	{
	if(!i)a=1.;
	else
		a=-i;
	b=x;
	podd= a*podd+b*peven;
	qodd= a*qodd+b*qeven;
	a= -(1. + i);
	b=1.;
	peven=a*peven+b*podd;
	qeven= a*qeven+b*qodd;
	zold=z;
	z=peven/qeven;
	itktg=i;
	if( (i>0) && abs(z-zold) < tol)
		{printf(" leaving eicf Ei/exp(x)=%e\n",z);
		return z*exp(x);
		}
	/* re norm*/
	if( abs(qodd)>10.)
		{
		y=1./qodd;
		peven*=y;
		podd*=y;
		qeven*=y;
		qodd=1.;
		}
	}
printf(" Ei not converged\n");
return(z*exp(x));/* return best guess*/
}
/*
double eia(x) double x;
{
double y,e1(),exp();
y=1./x;
if(y>.4)return(0.);
return e1(x)*exp(2.*x)*(1.+8.808*y-32.06*y*y+185.633*y*y*y)
	/(1.+6.5418*y-37.207*y*y+130.691*y*y*y);
}
*/

double e(x,n) double x;int n;
{
int k,ncross=2;
/* choose ncross based on timing*/
double y,ee;
if(!n) return exp(-x)/x;
if(n==1)return e1(x);
if(x>1. && n>ncross)return en(x,(double) n);
y=exp(-x);
ee=e1(x);
for(k=1;k<n;k++)
	{
	ee= (y-x*ee)/k;
	}
return ee;
}


double E1(x) double x;
{double Ei();
if(x<0.)return -Ei(-x);
if(x==0.)return -errorcode;
if(x<=1.)return e1s(x);
/* x>1.*/
return en(x,1.);
}



#define itmax 40

static double fcf(x,n,a,b) double x,a[],b[];int n;
{double value;int i;
value=1.;
for(i=1;i<=n;i++)
	{
	value= b[n-i]/(a[n-i]+x+value);
	}
return value;
}

static double a24[9]={-3.00000000320981265753,-5.00006640413131002475,
-7.06810977895039358836,-15.2856623636929636839,-7.63147701620253630855,
-27.9798538624305389340,-18.1949664929868906455,-223.12767077763240955,
175.33880126546597239},	b24[9]={
1.99999999999048104167,-2.99999894040324959612,-7.99243595776339741065,
-12.0187763547154743238,70.4831847180424675988,117.179220502086455287,
137.7903902365747998793,3.97277109100414518365,39784.597716741472084};

static double a12[8]={-2.073093182550626,66.81633520851786,
-61.8811145837267426,13.615334713984657547,-32.460770029937463678,
-8.7163735593963354058,
-17.8083116036997799966,-4.088838379362196702},b12[8]={
 .99905385353462753131,.99655264231071911439,4267.5815993503950849,
 11.526405585745173857,688.27350646689188421,27.549295584621895224,
 77.697693140151074176,50.693912820579735193};

static double a6[8]={5.731167057445080,4.1810242256285662231,5.886582407532811,
-19.41329675144307,7.8947220929445722122,23.273023383903914,
-36.778311347831145794,-2.4694098344836126512},b6[8]={
1.14625253249101619143,-199.1496002312351636,341.36521252437553905,
52.316556873455861379,317.27948925436932786,-8.38767084118964070656,
965.40521742928030312,2.6398300731802459334};


#define tol 1.e-6

static double Eict(x)double x;
{
double exp(),y,z;
if(x>24.)
	z=1.+(1.000000000000004855+fcf(x,9,a24,b24))/x;
else if(x<=24. && x>12.)
	z=(1.0000051738331117+fcf(x,8,a12,b12));
else if(x<=12. && x>=6.)
	{
	z=(9.9895766651165517e-01+fcf(x,8,a6,b6));
	}
	y=exp(x)/x;
/*printf(" exp/x=%e other %e ",y,z);*/
	return y*z;
}

double scale;

/* approximate accuracy: 
order	digits
3	5				
4	7
6	12 
9	20
*/
/*#define order 3*/
/*#define order 4*/
/*#define order 6*/
#define order 9

static double Eis(double x)
	{
	double x0= .372507410781366634461991866580;

/*	double p[4]={-2.2409438e3*.5,-1.832677e2,-9.0053723e1,-5.6796891};
	double q[4]={-8.5609557e2*.5, 4.097087e2,-8.3056797e1,6.75};
*/
/*	double p[5]={3.360281229e4*.5,1.521253684e3,1.540472598e3,
7.554503527e1,6.214595703};
	double q[5]={.5*1.337025876e4,-6.737936555e3,1.583666607e3,
-1.921193866e2,1.01250e1};
*/
/*	double p[7]={1.45747321743484e7*.5,4.71761893143542e4,7.84689687508871e5,
1.35036011184167e4,7.35420513947427e3,1.80339138516268e2,7.75745083328328};
	double q[7]={.5* 6.06704887208133e6,-3.22441109528198e6,8.67868680839173e5,
-1.39638471189123e5,1.39851277161255e4,-8.27775609141545e2,2.278125e1};
  */
double p[10]={-4.1658081333604994241879e11*.5,1.2177698136199594677580e10,
-2.5301823984599019348858e10,3.1984354235237738511048e8,
-3.53778096944311334848e8,-3.139866086424726586205e5,
-1.4299841572091610380064e6,-1.42870725001970057773776e4,
-1.28312206592620678144e3,-1.296370260247483002859e1};
double q[10]={-1.7934749837151009723371e11*.5,9.8900934262481749439886e10,
-2.8986272696554495342658e10,5.4229617984472955011862e9,
-7.0108568774215954065376e8,6.4698830956576428587653e7,
-4.2648434812177161405483e6,1.9418469440759880361415e5,
-5.5648470543369082846819e3,7.688671875e1};

	double h,z,y,n,no,noo,d,dold,doo;  int i;
	z=x/6.;y =2.*z-1.;z=2.*y;
	noo=doo=0.; no=p[order];dold=q[order];
	for(i= (order-1);i>0;i--)
		{
		n= z*no-noo+p[i];
		d= z*dold-doo+q[i];
		doo=dold;dold=d;
		noo=no;no=n;
		}
	n=y*n-noo+p[0];
	d=y*d-doo+q[0];
	h=x-x0;
	return log(1.+h/x0)+(h)*n/d;
	}

double eias(x) double x;
{double term,oterm,sum;
int i;
sum=0; oterm=term=1.;
for(i=1;i<itmax;i++)
	{
	sum+=term;
	oterm=term;
	term*= i/x;
	if(abs(term)>abs(oterm) || abs(term)<tol*abs(sum))break;
	}
/*printf(" sum=%e\n",sum);*/
return exp(x)/x*sum;
}

double Ei(x) double x;
{
int i;
double sum,term,tolr=1.e-7,xcross=25.;
/*"exact"*/
if(x<=0.)return errorcode;
if(x>xcross)return eias(x);
if(x>1.)return eicf(x);
sum= log(x)+ Egamma;
term=x;
for (i=1;i<itmax;i++)
	{
	sum+=term/(i);
	term*=x/(i+1);
	if( abs(term)<tol || abs(term)<abs(sum)*tolr)return sum;
	itktg=i;
	}
printf(" too many iterations\n");
return errorcode;
}

double ei(x) double x;
{
/* next 4 stmts for Rational approx. remove otherwise*/
if(x<=0.)return -e1(-x);
if(x>20.)return eias(x);
if(x>6.)return Eict(x);
return Eis(x);
}

double li(x) double x;
{
double ei(),log();
if(x<1.)return errorcode;
return ei(log(x));
}

double alpha(x,n) double x; int n;
{
return pow(x, -(double)(n+1))*
	gamma((double)(n+1))*(1.-incgam((double)(n+1),x));

/* first find alpha n=0, then use recurrence:*/
/*currently disabled, might be useful small n,x*/
/*{
z=exp(-x);y=z/x;
if(!n)return y;

for (k=1;k<=n;k++)
	{
	y=(z+k*y)/x;
	}
return y;
}*/
}

double beta(x,n) double x; int n;
{
double realn,exp(),eps=1.e-5,z,y,sinh,sign,bet,log();
int k;
if(x==0.)return 0.;
if(!n)
	{ /* beta sub-0(x)=2sinh(x)/x no divide by zero */
	z=exp(-x);y=1./z;return (y-z)/x;
	}
if(x<0. || n< 1)return errorcode;
realn=n;/* n==0 handled above so no log problems*/
if( x<(realn*.368+.821+.184*log(realn)))
	{/*backward recursion*/
printf(" backward\n");
	k=10+n;
	z=exp(-x);y=1./z;
	bet=0.;
	sign= (k%2)?-1.:1.;
	for( ;k>=n;k--)
		{
		bet= (x*bet+z+sign*y)/(k+1);
		sign=-sign;
		}
	return bet;
	}
/* forward recursion*/
sign=-1.;
if(x<eps)
	{
	x=eps;
	z= exp(-x);y=1./z;
	sinh=.5*(y-z);
	bet=2. *(1.+x*x*.1666666);
	}
else
	{
	z=exp(-x); y=1./z;
	sinh=.5*(y-z);
	bet=sinh*2./x;
	}
if(!n)return bet;
for (k=1;k<=n;k++)
	{
	bet=(y*sign-z+k*bet)/x;
	sign*=-1.;
	}
return bet;
}


double f01(x) double x;
{double log(),y,z;
y=(Egamma+log(x));      z=x*x;
return .5*y*y-pi*pi/24.-z*
(.125+z*(-2.60416632e-3+z*(3.857955e-5+z*(-3.87035e-7+z*2.61455e-9))));
}
double g01(x) double x;
{double log(),y,z;
y=(Egamma+log(x));      z=x*x;
return -.5*pi*y+x*
(.9999999992+z*(-.0185185136+z*(3.3332344e-4+z*(-4.0422785e-6+z*3.201246e-8))));
}

double c01(x) double x;
{double sin(),cos(),f01(),g01(),z;
if(x<0.)return errorcode;
if(x<=2.)return f01(x)*cos(x)+g01(x)*sin(x);
z=x*x;
return -(z*(20.96173922+z*(7.150169966+z*(.4371420242+z*.005415719186)))-3.85642854)
 /((100.+z*(70.34218899+z*(11.16783932+z*(.496599058+z*.005415884237))))*z);
}

double s01(x) double x;
{double sin(),cos(),f01(),g01(),z;
if(x<0.)return errorcode;
if(x<=2.)return f01(x)*sin(x)-g01(x)*cos(x);
z=x*x;
return -(.5423489064+z*(43.88413692+z*(14.246573437+z*(.8460999768+z*.009972860283))))
 /((100.+z*(60.93657+z*(8.503549177+z*(.3369899234+z*.003314700857))))*x*z);
}


double g2(x,y)double x,y;
{
double sqrt(),log(),atan(),xs,ys,p,q,r,sum,a,sa,f(),g();
double aa[3]={17.420076,.501312744,3.43966581};
double A[3]={.05299936,.422384803,.241865419};
int i;
xs=x*x;ys=y*y;
p=f(x*.5);q=g(x*.5);                 r=log(2.*y/x);
sum=0.;
for(i=0;i<3;i++)
	{
	a=aa[i];sa=sqrt(a);
	sum+=A[i]/(a+xs)*(.5*x*log((ys+a)/(a+.25*xs))-sa*atan((2.*y-x)*sa/(x*y+2.*a)));}
return -(p*q)+f(x)*r-.282750417*(r/x)-sum;

}
double f2(x,y)double x,y;
{
double sqrt(),log(),atan(),xs,ys,p,q,r,sum,a,sa,f(),g();
double aa[3]={21.850456,.770345382,4.55715659};
double A[3]={.163725227,.34115997,.428765629};
int i;
xs=x*x;ys=y*y;
p=f(x*.5);q=g(x*.5);                 r=log(2.*y/x);
sum=0.;
for(i=0;i<3;i++)
	{
	a=aa[i];sa=sqrt(a);
	sum+=A[i]/(a+xs)*(.5*log((ys+a)/(a+.25*xs))+x/sa*atan((2.*y-x)*sa/(x*y+2.*a)));}
return .5*(q*q-p*p)-g(x)*r+ .066349174*(2.-x/y+r)/xs+sum;
}

double c2(x,y) double x,y;
{double f(),g(),f2();
if(x<2.)return errorcode;
if(1.<=y && y<=x) return f2(x,y);
if(y<1.)return g(y)*g(x-y)-f(y)*f(x-y)-f2(x,x-y);
return errorcode;
}

double s2(x,y) double x,y;
{double f(),g(),g2();
if(x<2.)return errorcode;
if(1.<=y && y<=x) return g2(x,y);
if(y<1.)return -g(y)*f(x-y)-f(y)*g(x-y)-g2(x,x-y);
return errorcode;
}


struct complex ceiscaled;

/* continued fraction approx. from CALGO 14 */
cei(x,y,k,toler,u,v,n) double x,y,k,toler,*u,*v;int *n;
{
struct complex z,ans,e,temp; double d;
z.x=x;z.y=y;
d=cabs(z);
if(d<1.e-20)
	{*n=0;
	if(k>1){*u=1./(k-1.);*v=0.;}
	else {*u=-errorcode;*v=0.;}
	return 0;
	}
ceii(x,y,k,toler,u,v,n);
/*printf(" %e %e scaled zexp(z)E\n",*u,*v);*/
ceiscaled.x=*u;
ceiscaled.y=*v;
/* scale by dividing by z*exp(z)*/
CMPLX(ans,*u,*v);CDIV(temp,ans,z);cexp(&z,&e);
CDIV(ans,temp,e);
*u=ans.x;*v=ans.y;
return 0;
}


ceii(x,y,k,toler,u,v,n) double x,y,k,toler,*u,*v; int *n;
{/* finds complex 'exponential integral' u+iv= z^kexp(z) integral z to infinity
(z=x+iy) of exp-t dt/t^k . do not use |z|<.05 as convergence is slow*/
 double M,K,t1,t2,t3,a,b,c,d,g,h,e; int m;
 struct complex  cz,term,sum,power,offset;double fac;
 cz.x=x;cz.y=y;
 fac=cabs(cz);
 if(fac<1.e-20)
	{ if(k>1) {*u=0.;*v=0.;/* En=1./(k-1),but scaled x*exp(x)=0*/
		return;}
	  else {*u=-errorcode;*v=0.;return;}
	}
 if(fac<=.05 &&!(y==0.&&x<0))
	{/* use series for small |z| unless arg(z)=pi*/
	if(k==1)
		{CMPLX(offset,-Egamma,0.);
		clog(&cz,&term);CSUB(offset,offset,term);
		}
	else
		{
		clog(&cz,&term);CMPLX(offset,digamma((double)k),0.);
		CSUB(offset,offset,term);
		CTREAL(offset,offset,1./gamma((double)k));
		CTREAL(term,cz,-1.);clog(&term,&sum);
		CTREAL(sum,sum,(double)(k-1));cexp(&sum,&term);
		CMULT(sum,offset,term);CLET(offset,sum);
		}
	CMPLX(sum,0.,0.);
	CTREAL(cz,cz,-1.);
	fac=1.; CMPLX(power,1.,0.);
	for(m=0;m<100;m++)
		{
		if(m) fac*=m;
		if(m==(k-1))continue;
		CTREAL(term,power,1./((m-k+1)*fac));
		CADD(sum,sum,term);
		if(cabs(term)<cabs(sum)*toler)break;
		CMULT(term,power,cz);
		CLET(power,term);
		*n=m;
		}
	CSUB(offset,offset,sum);/* now multiplby by exp(x)*x */
	CMULT(sum,offset,cz);CTREAL(sum,sum,-1.);
	CTREAL(offset,cz,-1.);cexp(&offset,&term);
	CMULT(offset,term,sum);*u=offset.x;*v=offset.y;return;
	}
 e=tol*tol;
 *u=c=a=1.;*v=d=b=0.;
 *n=1; K=k-1;
 do
 {
 g=*u;h=*v;(*n)++;
/*printf(" n is %d\n",*n);*/
 m=(*n)>>1;
 M=((m<<1)==(*n))? m+K:m;
 t1=x+M*c;t2=y+M*d;
 if(t1==0. && t2==0.){fprintf(stderr," cei: t1=t2=0\n");return;}
 t3=1./(t1*t1+t2*t2);
 c=(x*t1+y*t2)*t3;
 d=(y*t1-x*t2)*t3;
 t1=c-1;t2=a;
 a=a*t1-d*b;b=d*t2+t1*b;
 *u=g+a;*v=h+b;
 }while( (a*a+b*b)/(*u* *u+*v* *v)>e);
 return;
}

cexpint(z,n,toler,ans,iter) struct complex *z,*ans;double toler,n;int *iter;
{ cei(z->x,z->y,n,toler,&(ans->x),&(ans->y),iter);
return 0;
}
