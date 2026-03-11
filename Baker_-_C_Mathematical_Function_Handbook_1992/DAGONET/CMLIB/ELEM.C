/* elementary transcendental functions
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

ln	natural logarithm
expon	exponential
sine 	sine functions
arc_tan	arc tangent
arc_sin	arc sine
	hyperbolic functions:
hyper_sin
hyper_cos
hyper_tan
	arc hyperbolic functions:
arc_hyper_sin
arc_hyper_cos
arc_hyper_tan

cosine
arc_cosine
tangent 	tangent= sin/cos
tangnt		tangent by continued fraction method
arc_tangent

cosine
*/
#include "cmlib.h"
#include "protom.h"
#define tol 1.e-9
#define twopi 2.*pi
#define halfpi .5*pi

double ln(double xin)
	{double sum,x,term,/*old,*/y,ys,factor,power,offset;int i;
	double e=2.718281828459045;
	/*positive x only, please!*/
	x=xin;
	if(x<=0.)return errorcode;
	if(x==1.)return 0.;
	if(x<.5)return -ln(1./x);
	offset=0.;
	power=factor=term=1.;
	sum= 1.;i=3;
	/* unsophisticated range reduction*/
	while( x>e)
		{  x/=e; offset+=1.;
		 /*printf("offset now %le x=%le\n",offset,x);*/
		}
	/* loop won't work well for large |x|, so need the above */
	y= (x-1.)/(x+1.);
	ys=y*y;
	infinite_loop
		{
		/*old=term;*/
		factor=1./((double)i);
		power*=ys;
		term=factor*power;
		sum+=term;
		if( abs(term)< tol)break;
		i+=2;
		}
	return 2.*y*sum+offset;
	}

double expon(double x)
	{
	double sum,term/*,old,y*/,factor,power;int i;
	if(x<-65.)return 0.;
	if(x>65.) return -errorcode;
	if(x>=1.)return 1./expon(-x);
	power=factor=term=1.;
	sum= 1.;i=1;
	infinite_loop
		{
		/*old=term;*/
		factor /= ((double)i);
		power*=x;
		term=factor*power;
		sum+=term;
		if( abs(term)< tol)break;
		i++;
		}
	return sum;
	}

double sine(double xi)
	{
	double sum,term,/*old,*/x,y,factor,power,z;/*int i;*/
	x=xi;
	while(x > twopi) x-=twopi;
	while(x < -twopi) x+=twopi;
	if(x==0.)return 0.;
	if(x==twopi)return 0.;
	power=factor=term=1.;
	sum= 1.;z=2.;y=x*x;
	infinite_loop
		{
		/*old=term;*/
		factor /=(-z*(z+1.));
		power*=y;
		term=factor*power;
		sum+=term;
		if( abs(term)< tol)break;
		z+=2.;
		}
	return sum*x;
	}


double arc_tan(double xi)
	{
	double sum,term,old,x,y,power/*,sign,factor*/,z;int i/*,big*/;
	x=xi;
	if(x==0.)return 0.;
	/* convergence problems if y=x=1.*/
/*  series  4.4.42 first
	if(x==1.)return pi*.25;
	if(x==-1.)return -pi*.25;
	power=factor=term=1.;
	sum= 1.;z=3.;
	if(abs(x)<1.)big=0;
	else {big=1;x=1./x;}
	y=x*x;          sign=-1.;
	infinite_loop
		{
		factor =(sign/z);
		power*=y;
		term=factor*power;
		sum+=term;
		if( abs(term)< tol)break;
		z+=2.;                sign=-sign;
		}
	sum=big?pi*.5-sum*x:sum*x;
	if(x>=0.)return sum;
	return -sum;
	*/
	if(abs(x)>1.)
		{ sum=arc_tan(1./x);
		if(x>0.)return pi*.5-sum;
		else return -pi*.5-sum;
		}
	/* series 4.42.last*/
	y=x*x;old=1./(1.+y);z=y*old;
	sum=1.;i=4;power=z*2./3.;
	infinite_loop
		{
		term=power;
		sum+=term;
		if( abs(term)< tol)break;
		power*=(z*i/(i+1.));
		i+=2;
		}
	return x*old*sum;
	}

double arc_sine(double xi)
	{
	double sum,term/*,old*/,x,y,factor,power,z,q;/*int i;*/
	x=xi;
	if(x==0.)return 0.;
	if(x==1.)return halfpi;
	if(x==-1.)return -halfpi;
	if(abs(x)>1.)return errorcode;
	if(abs(x)>.7071)
		{
		sum=arc_tan(x/square_rt(1.-x*x));
		return sum;
		}
/*		return halfpi-arc_sine(square_rt(1.-x*x));*/
	power=factor=term=1.;
	sum= 1.;z=2.;y=x*x;
	while(1)
		{
		/*old=term;*/
		q=z-1.;
		factor *=    q*q/(z*(z+1.));
		power*=y;
		term=factor*power;
		sum+=term;
		if( abs(term)< tol)break;
		z+=2.;
		}
	return sum*x;
	}

double hyper_sin(double x)
	{double ex;ex=expon(x); return (ex-1./ex)*.5;}
double hyper_cos(double x)
	{double ex;ex=expon(x); return (ex+1./ex)*.5;}
double hyper_tan(double x)
	{double ex,i;ex=expon(x);
	if(ex!=0.)i=1./ex;
	else return 1.;
	return (ex-i)/(ex+i);}
double arc_hyper_sin(double x)
	{ return ln(x+square_rt(x*x+1.));}
double arc_hyper_cos(double x)
	{
	if(x<1.)return errorcode;
	if(x==1.)return 0.;
	return ln(x+square_rt(x*x-1.));}
double arc_hyper_tan(double x)
	{
	if(abs(x)>1.)return errorcode;
	if(x==1.)return -errorcode;
	if(x==-1.)return errorcode;
	return .5*ln((1.+x)/(1.-x));}

double cosine(double x) { return sine(pi*.5+x);}
double arc_cosine(double x) {return pi*.5-arc_sine(x);}

double tangent(double x)
	{
	double d;d=cosine(x);
	if(d==0.)return errorcode;
	return sine(x)/d;
	}

double tangnt(double x)
	{/*continued fraction method*/
	double z,ao=1.,ae=0.,bo=0.,be=1.,c,a,b,old,new,norm;int first=1;
	while(x>pi)x-=pi;
	while(x<-pi)x+=pi;
	if(abs(x)==halfpi)return errorcode;
	if(x==0.)return 0.;
	/* tan= z/1- z^2/3- z^/5-*/
	z=x*x;
	c=1.;
	old=0.;
	norm=1.;
	while(1)
		{
		a=first? x:-z;
		b=c;
		ao= b*ae+a*ao;
		bo= b*be+a*bo;
		c+=2.;
		a=-z;b=c;
		ae= b*ao+a*ae;
		be= b*bo+a*be;
		/*printf(" ae be %le %le\n",ae,be);*/
		if( abs(ae)>10.)
			{norm=1./(ae);
			ae*=norm;be*=norm;
			ao*=norm;bo*=norm;
			}
		c+=2.;
		first=0;
		new= ae/be;
		if(abs(new-old)<tol)break;
		old=new;
		}
	return new;
	}
/*continued fraction method slower than tangnt*/
/*
double Tangnt(double x)
	{
	double z,ao=1.,ae=0.,bo=0.,be=1.,c,a,b,old,new;
	while(x>pi)x-=pi;
	while(x<-pi)x+=pi;
	if(abs(x)==halfpi)return errorcode;
	if(x==0.)return 0.;
	z=x*x;
	c=1.;
	old=0.;
	while(1)
		{
		a=-z;
		b=c;
		ao= b*ae+a*ao;
		bo= b*be+a*bo;
		c+=2.;
		a=-z;b=c;
		ae= b*ao+a*ae;
		be= b*bo+a*be;
		c+=2.;
		new= ae/be;
		if(abs(new-old)<tol)break;
		old=new;
		}
	return -new/x;
	}
*/
double arc_tangent(double y, double x)
	{double at;
	/* returns angle between -pi and pi.  "branch cut" x<0 */
	if(x==0.)
		{if (y==0.)return 0.;/* or errorcode-this is choice of atan2*/
		if(y>0.)return halfpi;
		else return  -halfpi;
		}
	at=arc_tan( y/x);/* sign of at is sign of y/x -pi/2<at<pi/2*/
	/* positive at*/
	if(x>0. )return at;/* quadrant I or 4 -pi/2 to pi/2 */
	if(y<0.)return at-pi;/* or pi+at>pi quadrand III*/
	/* negative at*/
	/* else quadrant II*/
	return at+pi;
	}
