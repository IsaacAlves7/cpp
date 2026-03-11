/*  Clebsh-Gordon, Wigner, related coefficients
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

triangle	tests triangle inequality
DELTA		delta function (inline)
isint,notint	is the (type double) an integer?
m1e		-1^x

wigner()	wigner coef. efficient
		(minimal calls to gamma function)
Clebsh-Gordon	Clebsh-Gordon coef. calls wigner
wigner3j	3j symbol calls wigner.
wigner6j	wigner 6j calls racah.
Wigner6j	wigner 6j calls racah. more natural order of arg.
racah		W of Racah
V		V of Racah
X		X of Racah =Wigner9j
CG		Clebsh-Gordon calls V
Wigner9j	9j symbol uses Wigner6j

based upon: Racach Phys Rev 62 438 (1942) (W,V,X)
Angular Momentum, Brink & Satchler, Oxford U Press 1968 9j
Bethe & Jackwicz, Intermediate Quant. Mech & Abramowitz&Stegun, wigner
Mizushima, Theor. Phys. Wiley,1972. 9j,special values
 */
#include "cmlib.h"
#include "protom.h"
#include <stdio.h>

#define ab(x) ((x)>0. ? (x):0.)
#define min(a,b) ((a)<(b)? (a):(b))
#define max(a,b) ((a)>(b)? (a):(b))
#define DELTA(a,b) ((a)==(b)?1.:0.)
#define tol  3.e-7
#define gamret 1.e60

int notint(x) double x;
{ double y;
y=x-((int)x);
return y==0.?0:1;
}
/* check if multiple of 1/2 ideally,use to chk all arguments*/
int nothint(x) double x;
{ double y;
y=2.*x;
y=y-((int)y);
return abs(y)<tol?0:1;
}

double gam(x) double x;
{
int y;
y=x; if(!y)return gamret;
if( abs(x-(double)y)< tol && y<0)return gamret;
return gamma(x);
}

int triangle(j1,j2,j) double j1,j2,j;
{ if( abs(j1-j2)<=j && j<=(j1+j2))return 1;
return 0;
}

double m1e(x)double x;
{int j;
j=x;
if(x!=(double)j){fprintf(stderr," -1 to non-integer=%f\n",x);return 0.;}
return (j%2)?-1.:1.;
}

/* Wigner or Clebsh-Gordon coef.*/
double wigner(j,j1,j2,m,m1,m2) double j,j1,j2,m,m1,m2;
{double coef,arg1,arg2,arg3,arg4,arg5,sum,term;
int k,kmin,kmax;
if(abs(m-m1-m2)>1.e-8) return 0.;
if(j2==0. && m2==0.) return  (j==j1 && m==m1)? 1.: 0.;
if( abs(m1)>j1 || abs(m2)>j2 || abs(m)>j ||
 j1+j2-j<0. || j1-j2+j <0. || j+j2-j1<0.) return 0.;
coef= sqrt((2.*j+1.)*gam(j1+m1+1.)*gam(j1-m1+1.)*gam(j2+m2+1.)
 *gam(j2-m2+1.)*gam(j+m+1.)*gam(j-m+1.)
 *gam(j1+j2-j+1.)*gam(j+j1-j2+1.)*gam(j+j2-j1+1.)/gam(j+j1+j2+2.));

 kmin=max(0,max(j2-j-m1,j1-j+m2));
 arg1=gam(j1+j2-j+1.-kmin);
 arg2=gam(j1-m1+1.-kmin);
 arg3=gam(j2+m2+1.-kmin);
 arg4=gam(j-j2+m1+1.+kmin);
 arg5=gam(j-j1-m2+1.+kmin);
 if(arg1==errorcode ||arg2==errorcode||arg3==errorcode||arg4==errorcode
 || arg5==errorcode) sum=0.;
 else sum= m1e((double)kmin)/(arg1*arg2*arg3*arg4*arg5);
/*printf(" coef=%e\n",coef);
printf(" k=0 sum=%e\n",sum);*/
if(sum!=0.)term=sum;
else {term=1.;
	if(arg1!=errorcode) term*=arg1;
	if(arg2!=errorcode) term*=arg2;
	if(arg3!=errorcode) term*=arg3;
	if(arg4!=errorcode) term*=arg4;
	if(arg5!=errorcode) term*=arg5;
/*printf(" arg %e %e %e %e %e\n",arg1,arg2,arg3,arg4,arg5);*/
	term=1./term;
	}
/*printf(" term=%e\n",term);*/
kmax=min(j2+m2+1,min(j1+j2-j+1, j1+1-m1));
/*printf("kmin, kmax %d %d\n",kmin,kmax);*/
kmin++;
for(k=kmin;k<kmax;k++)
{
	arg1=(j-j2+m1+k);
	arg2=(j-j1-m2+k);
	if(arg1<=0.)arg1=1.;
	if(arg2<=0.)arg2=1.;
	term=-term*(j1+j2-j-k+1)*(j1-m1-k+1)*(j2+m2-k+1)/(arg1*arg2*k);
/*printf(" arg1,2=%e %e term=%e sum=%e \n",arg1,arg2,term,sum);*/

	sum+=term;
/*printf(" sum, term now %e %e\n",sum,term);*/
	if(abs(term)<=abs(sum)*tol)break;
	}
return sum*coef;
}

double ClebshGordon(j1,m1,j2,m2,j,m)double j,m,j1,m1,j2,m2;
{return wigner(j,j1,j2,m,m1,m2);}


/*wigner 3j coefficients*/
double wigner3j(j,j1,j2,m,m1,m2) double j,j1,j2,m,m1,m2;
{
if((m1+m2+m)!=0. || !triangle(j1,j2,j))return 0.;
return wigner(j,j1,j2,-m,m1,m2)/sqrt(2.*j+1.)*m1e(j1-j2-m);
}

double delta(a,b,c)double a,b,c;
{double value;
value= sqrt(ab(gam(a+b-c+1.)*gam(b+c-a+1.)*gam(a+c-b+1.)
	/gam(a+b+c+2.)));
return value;
}

double gammn(x) double x;
{if(x>0.)return gamma(x);return 0.;}

/* Racah W(abcd;ef)*/
double racah(a,b,c,d,e,f) double a,b,c,d,e,f;
{/* arguments:j1 j2 j3 j j12 j23 */
double coef,sum,sign;int k;
double z,zmin,zmax,q;
/* See Racah Phys. Rev 62 eqn 39 p.444 for triad sum condition*/

if(nothint(a)||nothint(b)||nothint(c)||nothint(d)|nothint(e)||nothint(f))
	{fprintf(stderr," Racah W:argument not n/2 \n");return 0.;}
if(notint(a+b+e) ||notint(c+d+e)||notint(a+c+f)||notint(b+d+f))
	{fprintf(stderr," Racah W: nonintegral triaad sum\n");return 0.;}
if(!triangle(a,b,e) ||!triangle(c,d,e)||!triangle(b,d,f)||!triangle(a,c,f))
	return 0.;
coef=delta(a,b,e)*delta(a,c,f)*delta(b,d,f)*delta(c,d,e);
/*printf(" coef %e\n",coef);*/
if(coef==0.)return 0.;
sum=0.;
zmax=min( a+b+c+d+1.,min(c+d-e,min(a+b-e,min(a+c-f,(b+d-f)))))+1.;
/*zmin=  -max(e+f-a-d,e+f-b-c);*/
zmin=0.;
k=zmin;
sign= k%2?-1.:1.;
if( zmin==(double)k && zmin<0.){zmin=0.;sign=1.;}
/*printf(" zmin, zmax %e %e sign %e\n",zmin,zmax,sign);*/
for(z=zmin;z<=zmax;z+=1.)
	{sum+=q=sign*gammn(a+b+c+d+2.-z)/
		(gam(z+1.)*gam(e+f-a-d+z+1.)*gam(e+f-b-c+z+1.)
		*gam(a+b-e-z+1.)
		*gam(c+d-e-z+1.)*gam(a+c-f-z+1.)*gam(b+d-f-z+1.));
	if(abs(sum)<1.e-20 && zmin<0.)sum=0.;
	else sign=-sign;
	/*printf(" sum now %e z was %e term %e\n",sum,z,q);*/
	}

/*printf(" coef, sum %le %le\n",coef,sum);*/
return coef*sum;
}


/*  form of wigner 6-j symbol {  a b e  }
				 c d f          */
double wigner6j(a,b,c,d,e,f) double a,b,c,d,e,f;
{return (((int)(a+b+c+d))%2?-1.:1.)*racah(a,b,d,c,e,f);}


/*  form of Wigner 6-j symbol {  a b c  }
				 d e f
*/
double Wigner6j(a,b,c,d,e,f) double a,b,c,d,e,f;
{return (((int)(a+b+e+d))%2?-1.:1.)*racah(a,b,e,d,c,f);}

/* form of Wigner 9-j: { a b c }
			 d e f
			 g h i
*/
double Wigner9j(a,b,c,d,e,f,g,h,i)double a,b,c,d,e,f,g,h,i;
{double sum,k,coef,kmin,kmax,kinc,q;
sum=0.;
kmax= abs(a)+abs(b)+abs(c)+abs(d)+abs(e)+abs(f)+abs(g)+abs(h)+abs(i);
kmin=0.;kinc=1.;
for(k=kmin;k<=kmax;k+=kinc)/*increment k by one? kmin =0,1, or smth else?*/
	{
	coef=1.;if((2.*k-((int)(2.*k)))!=0.)coef=-1.;
	sum+=q=(2.*k+1.)*coef*Wigner6j(a,b,c,f,i,k)
	*Wigner6j(d,e,f,b,k,h)*Wigner6j(g,h,i,k,a,d);
	}
return sum;
}

/* X should be Wigner9j*/
double X(a,b,c,d,e,f,g,h,i)double a,b,c,d,e,f,g,h,i;
{double sum,k,coef,kmin,kmax,kinc,q;
sum=0.;
kmax= abs(a)+abs(b)+abs(c)+abs(d)+abs(e)+abs(f)+abs(g)+abs(h)+abs(i);
/*kmax*=3.;*/
kmin=0.;/* was -kmax*/
kinc=1.;
for(k=kmin;k<=kmax;k+=kinc)
	{
	if(!triangle(k,a,i)||!triangle(k,d,h)||!triangle(k,b,f))continue;
	sum+=q=(2.*k+1.)*racah(a,i,d,h,k,g)
	*racah(b,f,h,d,k,e)*racah(a,i,b,f,k,c);
	/*printf(" X k,q sum %le %le %le\n",k,q,sum);*/
	}
return sum;
}

/* V, CG based on Racah method in Bethe & Jacwitz*/
/* Racach V(a,b,c;A,B,C) */
double V(a,b,c,A,B,C) double a,b,c,A,B,C;
{double term,coef,sum,factor,p,q;int k;
if(notint(a+A)||notint(a-A)||notint(b+B)||notint(b-B)||
	notint(c+C)||notint(c-C)||
	notint(a+b-c)||notint(a+c-b)||notint(b+c-a))
	{fprintf(stderr," nontinegral sums in Racah V\n");return 0.;}
if((a+A)<0.||(a-A)<0.||(b+B)<0.||(b-B)<0.||
	(c+C)<0.||(c-C)<0.||
	(a+b-c)<0.||(a+c-b)<0.||(b+c-a)<0.){return 0.;}
/*if(abs(A)>a||abs(B)>b||abs(C)>c)return 0.;*/
coef= gam(a+b-c+1.)*gam(a-b+c+1.)*gam(b+c-a+1.)/gam(a+b+c+2.)
	*gam(a+A+1.)*gam(a-A+1.)*gam(b+B+1.)*gam(b-B+1.)
	*gam(c+C+1.)*gam(c-C+1.);
sum=0.;factor=1.;
for(k=0;k<1000;k++)
	{
	p=1.+k;q=1.-k;
	term=factor/(gam(p)*gam(a+b-c+q)*gam(a-A+q)
	*gam(b+B+q)*gam(c-b+A+p)*gam(c-a-B+p));
	sum+=term;factor=-factor;
	if(abs(term)< 1.e-6*abs(sum))break;
	}
return sqrt(coef)*sum;
}

double CG(j1,j2,j,m1,m2,m)double j1,j2,j,m1,m2,m;
{
if(m!=(m1+m2)) return 0.;
return sqrt(2.*j+1.)*V(j1,j2,j,m1,m2,m);
}
