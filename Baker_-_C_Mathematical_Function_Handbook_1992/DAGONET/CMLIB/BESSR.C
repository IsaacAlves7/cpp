/*  
Bessel functions: rational approximations
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

limit	error magnitudes:  |x|<=3  	|x|>3
j0									5e-8	 	1.6e-8
y0 							   1.4e-8 	1.6e-8
j1                         1.3e-8 	1.1e-8
y1                         1.1e-7 	1.1e-8

limit	error magnitudes:  |x|<=3.75  |x|>3.75
i0					            1.6e-7   1.9e-7
i1                         8.e-9    2.2e-7

limit	error magnitudes:   |x|<=2    |x|>2
k0                         1.e-8    1.9e-7
k1									8.e-9    2.2e-7
be,ke kelvin functions and derivatives, 0th order.
limit	error magnitudes:  |x|<=8.    |x|>8
ber                        1.e-9     3.e-7
bei                        1.e-9     3.e-7
ker                        1.e-8     1.e-7
kei                        3.e-9     1.e-7
ber'                       2.1e-8    3.e-7
bei'                       7.e-8     3.e-7
ker'                       8.e-8     2.e-7
kei'                       7.e-8     2.e-7

*/

#include "cmlib.h"
#include "protom.h"

double j0(x) double x;
{	
double t,cos(),sqrt(),f,y,bj;
if(x<3.)
	{
	y=(x*.33333333);y=y*y;
	bj= (((((.00021*y-.0039444)*y+.0444479)*y-.3163866)*y
        +1.2656208)*y-2.2499997)*y+1.;
	return(bj);
	}
y=3./x;	
if(x<1.e6)
	{
	f=(((((.00014476*y-.00072805)*y+.00137237)*y-.00009512)*y
     -.00552740)*y-.00000077)*y+.79788456;
	t=(((((.00013558*y-.00029333)*y-.00054125)*y+.00262573)*y
     -.00003954)*y-.04166397)*y-.78539816+x;
	bj=f*cos(t)/sqrt(x);
	return(bj);
	}
	return(0.);
	}

double y0(x) double x;
{
double f,t,y,log(),sin(),sqrt(),j0(),by;
if(x<3. && x>=0.)
	{
	y=(x*.333333333);y=y*y;
	by=(((((-.00024846*y+.00427916)*y-.04261214)*y+.25300117)*y
     -.74350384)*y+.60559366)*y+.36746691 
     +j0(x)*log(x*.5)*.636619772;
	return(by);
	}
if(x>=3.)
	{
	y=3./x;
	if(x> 1.e6)return(0.);
	f=(((((.00014476*y-.00072805)*y+.00137237)*y-.00009512)*y
       -.00552740)*y-.00000077)*y+.79788456;
	t=(((((.00013558*y-.00029333)*y-.00054125)*y+.00262573)*y
      -.00003954)*y-.04166397)*y-.78539816+x;
	return(f*sin(t)/sqrt(x));
	}
return(-1.);
}

double j1(x)double x;
{
double y,t,f,sqrt(),cos(),bj;
if(x<=3. && x>0.)
	{
	y=(x*.333333333);y=y*y;
	bj= ( (((((.00001109*y-.00031761)*y+.00443319)*y-.03954289)*y
      +.21093573)*y-.56249985)*y+.5)*x;
	return(bj);
	}
if(x>3.)
	{
	if(x>1.e6)return(0.);
	y=3./x;
	f=(((((-.00020033*y+.00113653)*y-.00249511)*y+.00017105)*y
      +.01659667)*y+.00000156)*y+.79788456;
	t=(((((-.00029166*y+.00079824)*y+.00074348)*y-.00637879)*y
      +.0000565)*y+.12499612)*y-2.35619449+x;
	return(f*cos(t)/sqrt(x));
	}
return(0.);
}

double y1(x) double x;
{
double y,j1(),log(),sin(),t,f,by;
if(x<3. && x>0.)
	{
	y=(x*.33333333);y=y*y;
	by=((((((.0027873*y-.0400976)*y+.3123951)*y-1.3164827)*y
     +2.1682709)*y+.2212091)*y-.6366198)/x
      +j1(x)*log(x*.5)*.636619772;
	return(by);
	}
y=3./x;
if(x>1.e6)return(0.);
	f=(((((-.00020033*y+.00113653)*y-.00249511)*y+.00017105)*y
     +.01659667)*y+.00000156)*y+.79788456;
	t=(((((-.00029166*y+.00079824)*y+.00074348)*y-.00637879)*y
     +.0000565)*y+.12499612)*y-2.35619449+x;
	return(f*sin(t)/sqrt(x));
}



double i0(x) double x;
{
double t,exp(),sqrt(),y;
if(x<3.75)
	{
	t=(x/3.75);t=t*t;
	y=((((((.0045813*t+.0360768)*t+.2659732)*t+1.2067492)
     *t+3.0899424)*t+3.5156229)*t+1.);
	return(y);
	}
t=3.75/x;
	y=((((((((.00392377*t-.01647633)*t+.02635537)*t
 	-.02057706)*t+.00916281)*t-.00157565)*t+.00225319)*t
    +.01328692)*t+.39894228)*exp(x)/sqrt(x);
    return(y);
}

double i1(x)double x;
{
double t,y,exp(),sqrt();
if(x<3.75)
	{
	t=x/3.75;t=t*t;
	y=((((((.00032411*t+.00301532)*t+.02658733)*t
      +.15084934)*t+.51498869)*t+.87890594)*t+.5)*x;
	return(y);
	}
t=3.75/x;
	y=((((((((-.00420059*t+.01787654)*t-.02895312)*t+
      .02282967)*t-.01031555)*t+.00163801)*t-.00362018)*t
     -.03988024)*t+.39894228)*exp(x)/sqrt(x);
     return(y);
}
	
double k0(x)double x;
{
double t,log(),i0(),y;
if(x<2.)
	{
	t=(.5*x);t=t*t;
	y=((((((.0000074*t+.0001075)*t+.00262698)*t+
     .0348859)*t+.23069756)*t+.4227842)*t-.57721566)
     -log(.5*x)*i0(x);
     return(y);
     }
 t=2./x;
	y=((((((.00053208*t-.0025154)*t+.00587872)*t
     -.01062446)*t+.02189568)*t-.07832358)*t+1.25331414)	
     *exp(-x)/sqrt(x);
     return(y);
}
double k1(x) double x;
{
	double y,t,exp(),log(),sqrt();
if(x< 2.)
	{t=(.5*x);t=t*t;
	y=((((((-.00004686*t-.00110404)*t-.01919402)*t
     -.18156897)*t-.67278579)*t+.15443144)*t+1.)/x
     +log(.5*x)*i1(x);
     return(y);
	}
t=2./x;
	y=((((((-.00068245*t+.00325614)*t-.0078353)*t
     +.01504268)*t-.03655620)*t+.23498619)*t+1.25331414)
      *exp(-x)/sqrt(x);
      return(y);
}


/* programs for computations of kelvin functions and their
	derivatives
*/

void ke(x,ker,kei,kerp,keip)double x,*ker,*kei,*kerp,*keip;
{
void be();
double y,z,al,sqrt(),log(),ber,bei,berp,beip;
struct complex CC,cex,cdum,f,phi,theta;
if(x<=8.)
	{y=x*x/64.;
	z=y*y;
	be(x,&ber,&bei,&berp,&beip);
	al=-log(.5*x);
	*ker=al* ber+.7853981634* bei
     +(((((((-.00002458*z+.00309699)*z-.19636347)*z+5.65539121)*z
     -60.60977451)*z+171.36272133)*z-59.05819744)*z-.57721566);
	*kei=al* bei-.7853981634* ber
     +(((((((.00029532*z-.02695875)*z+1.17509064)*z-21.30060904)*z
     +124.2356965)*z-142.91827687)*z+6.76454936)*y);
	*kerp=al* berp- ber/x+beip*.7853981634
     +x*((((((-.00001075*z+.00116137)*z-.06136358)*z+1.4138478)*z
     -11.36433272)*z+21.42034017)*z-3.69113734)*y;
	*keip=al* beip- bei/x-.7853981634* berp
     +x*((((((.00011997*z-.00926707)*z+.33049424)*z-4.65950823)*z
     +19.41182758)*z-13.39858846)*z+.21139217);
	return;
	}
/*else*/

	CMPLX(CC,-.7071067812,-.7071067812);
	CTREAL(CC,CC,x);
	CTHET(-x,&y,&z);
	CMPLX(theta,y,z);
	CADD(CC,theta,CC);
	cexp(&CC,&cex);
    y=1.253314137/sqrt(x);
    CTREAL(f,cex,y);
	*ker=f.x;
	*kei=f.y;
	 CPHI(-x,&y,&z);
	CMPLX(phi,y,z);/*cex now phi of AS*/
	CMULT(cdum,f,phi);
	*kerp=-cdum.x;
	*keip=-cdum.y;
return;
}

void be(X,BER,BEI,BERP,BEIP)double X, *BER, *BEI, *BERP, *BEIP;
{
struct complex CC,cex,cdum,g,theta,phi;
double Y,Z,sqrt(),FKER,FKEI,FKERP,FKEIP;
static double PII=.3183098862;
if(X<=8.)
	{
	Y=(X/8.);Y=Y*Y;
	Z=Y*Y	;
	*BER=((((((-.00000901*Z+.00122552)*Z-.08349609)*Z
      +2.64191397)*Z-32.36345652)*Z+113.77777774)*Z-64.)*Z+1.;
	*BEI=Y*((((((.00011346*Z-.01103667)*Z+.52185615)*Z
      -10.56765779)*Z+72.81777742)*Z-113.77777774)*Z+16.);
	*BERP=X*((((((-.00000394*Z+.00045957)*Z-.02609253)*Z
     +.66047849)*Z-6.06814810)*Z+14.22222222)*Z-4.)*Y;
	*BEIP=X*((((((.00004609*Z-.00379386)*Z+.14677204)*Z
     -2.31167514)*Z+11.37777772)*Z-10.66666666)*Z+.5);
	return;
    }
/*else*/

	CMPLX(CC,.7071067812,.7071067812);
	CTREAL(CC,CC,X);
	CTHET(X,&Y,&Z);
	CMPLX(theta,Y,Z);
	CADD(CC,theta,CC);
	cexp(&CC,&cex);
	Y=.3989422804/sqrt(X);
    CTREAL(g,cex,Y);
	ke(X,&FKER,&FKEI,&FKERP,&FKEIP);
	*BER=g.x-PII*FKEI;
	*BEI=g.y+PII*FKER;
	CPHI(X,&Y,&Z);
	CMPLX(phi ,Y,Z);
	CMULT(cdum,phi,g);
	*BERP= cdum.x-PII*FKEIP;
	*BEIP= cdum.y+PII*FKERP;

return;
}

CTHET(X,PARTR,PARTI)double X,*PARTI,*PARTR;
{
double Y;
	Y=8./X;
	*PARTI=((((.0000019*Y+.0000051)*Y*Y-.0000901)*Y-.0009765)
      *Y-.0110485)*Y-.3926991;
	*PARTR=((((.0000006*Y-.0000034)*Y-.0000252)*Y-.0000906)*Y*Y
	+.0110486)*Y;
return;
}
CPHI(X,PARTR,PARTI)double X,*PARTR,*PARTI;
{double Y;
	Y=8./X;
	*PARTI=(((((-.0000032*Y-.0000024)*Y+.0000338)*Y+.0002452)*Y
     +.0013811)*Y-.0000001)*Y+.7071068;
	*PARTR=((((((.0000016*Y+.0000117)*Y+.0000346)*Y+.0000005)*Y
	-.0013813)*Y -.0625001)*Y+.7071068);
     return;
}
