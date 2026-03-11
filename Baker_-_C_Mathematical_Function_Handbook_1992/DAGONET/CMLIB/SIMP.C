/* 
adaptive integration routine 
uses Simpson's rule adaptively
(from "C Tools for  Scientists and Engineers" by L. Baker)

CONTENTS:

adsimp
	adaptively computes integral of function

simp
	recursively applies simpson's rule

main
	test driver

funct
	function to be integrated

DEPENDENCIES:
NONE
*/


/* maxlev = max. allowed level
   levmax = max. level used
   level  = current level
   feval  = number of function evaluations used*/

/*#include "libc.h"
#include "math.h"
*/
#define abs(x) ((x)?(x):-(x))
#define max(a,b) ((a)>(b)?(a):(b))
int level,maxlev,levmax,feval;

double adsimp(a,b,eps,f) double a,b,eps,(*f)();
{
double simp(),aa,epss,absarea,ans,fa,fb,fm,range,est;
levmax=1;
feval=3;
level=1;
/*maxlev=6;*/
aa=a;epss=eps;
absarea=1.;
est=1.;
/*printf(" in adsimp %le %le\n",a,b);*/
range=b-a;
fa=(*f)(a);
fb=(*f)(b);
fm=(*f)(.5*(a+b)) *4.;
/* printf(" intergrand fa,fb,fm %e %e %e\n",fa,fb,fm);*/
ans=simp(aa,range,fa,fm,fb,absarea,est,epss,f);
/*printf(" in adsimp simp=%le\n",ans);*/
return(ans);
}

double simp(a,da,fa,fm,fb,area,est,eps,f) 
double a,da,fa,fm,fb,area,est,eps,(*f)();
{
double absarea;
double dx,x1,x2,est1,est2,est3,f1,f2,f3,f4,sum,epss,norm=.588;
absarea=area;
/*printf("simp %e %e %e %e %e \n %e %e %e %d\n"
,a,da,fa,fm,fb,absarea,est,eps,level);*/
dx=.333333333*da;
epss=eps*norm;
x1=a+dx;
x2=x1+dx;
f1=4.*(*f)(a+.5*dx);
f2=(*f)(x1);
f3=(*f)(x2);
f4=4.*(*f)(a+2.5*dx);
/*
f1=4.*funct(a+.5*dx);
f2=funct(x1);
f3=funct(x2);
f4=4.*funct(a+2.5*dx);
*/
feval+=4;
est1=(fa+f1+f2)*dx*.166666666666;
est2=(f2+fm+f3)*dx*.166666666666;
est3=(f3+f4+fb)*dx*.166666666666;
absarea=area-abs(est)+abs(est1)+abs(est2)+abs(est3);
sum=est1+est2+est3;
level++;
levmax= max(level,levmax);
if ( ((abs(est-sum)> eps*absarea)||(est==1.)) && level<maxlev)
	 sum= simp(a,dx,fa,f1,f2,absarea,est1,epss,f)
		+simp(x1,dx,f2,fm,f3,absarea,est2,epss,f)
		+simp(x2,dx,f3,f4,fb,absarea,est3,epss,f);

level--;
return(sum);
}
