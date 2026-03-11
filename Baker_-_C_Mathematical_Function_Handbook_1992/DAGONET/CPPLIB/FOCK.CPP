/*  Fock Functions

from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/

#include "complex.hpp"

#define pi 3.141592653589793238462643383279
#define errorcode -1.e60
#define ABS(x)  ((x)>0.?(x):-(x))

// -zeros of Ai
double zai[10]={2.33810741,4.08794944,5.52055983,6.78670809,7.94413,9.02265,
10.0402,11.0085,11.9300,12.8288};
// - zeros of Ai'
double zaip[10]={1.01879,3.24820,4.82009,6.16330,7.37217,8.48849,9.53545,
10.5277,11.4751,12.3848};
// values of Ai' at zeros of Ai
double aiz[10]={.70121,-.80311,.86520,-.91085,.94734,
-.97792281,1.00437012,-1.05776869,1.04872065,-1.06779386};
// values of Ai at zeros of Ai'
double aipz[10]={.53566,-.41902,.38041,-.35791,.34230,
-.33047623,.32102229,-.31318539,.30651729,-.30073083};

// w1= sqrt(pi)*(Bi+Ai) w2= sqrt(pi)*(Bi-Ai) v=-2j(w1-w2)=sqrt(pi)Ai w1=w2*

static double ga[28]={1.861,1.833,1.802,1.766,1.726,1.682,1.633,1.58,
1.523,1.463,1.399,1.344,1.266,1.197,1.113,1.059,.991,.925,.859,.797,.738,
.488,.315,.203,.13,.054,.022,.009};

static double gp[28]={15.43,10.07,5.78,2.47,0.07,358.49,357.63,357.42,
357.79,358.67,360.,1.71,3.75,6.07,8.62,11.36,14.24,17.24,20.32,
23.46,26.64,42.56,57.98,72.9,87.57,116.75,145.93,175.12};

static double fa[25]={2.160,1.992,1.829,1.672,1.521,1.377,1.241,1.112,
.992,.879,.776,.681,.594,.516,.446,.383,.327,.279,.236,.199,.167,
.0665,.025,.0091,.0033};

static double fp[25]={295.88,291.82,288.97,287.25,286.57,286.85,288.01,
289.98,292.68,296.04,300.,304.99,309.46,314.83,320.58,326.63,332.95,
339.49,346.22,353.09,.08,35.88,71.61,106.51,140.67};

int fockfg(double t, complex& f, complex& g)
{int i;double z,eps,mag,phase; complex factor,iu(0.,1.),sum;
if( t < -1.)
	{
	// asymptotic result for large negative t.
	// nothing better for range -1 to -2 available
	z=t*t*t;
	f= iu*2.*t* (-iu*z/3.).cexp()*(1.-iu*.25/z+.5/(z*z) );
	g= 2.* (-iu*z/3.).cexp()*(1.+iu*.25/z-1./(z*z) );
	mag= f.abs();phase= atan2( f.imaginary(),f.real());
	printf(" f mag=%le phase(rad)=%le\n",mag,phase);
	mag= g.abs();phase= atan2( g.imaginary(),g.real());
	printf(" g mag=%le phase(rad)=%le\n",mag,phase);
	return 1;
	}
eps=.00001;
if(t<=1.)
	{
	i= (t+1.+eps)/.1; z= (t-i*.1)+1.;
	}
else if(t<=3.)
	{i= 20+(t-1.+eps)/.5; z= t-(i-20)*.5-1.;
	}
else
	{i=24+(t-3.+eps);z=t-(i-24)-3.;
	}
//{printf("  z=%le t=%le i=%d\n",z,t,i);}
mag= ga[i]*(1.-z)+z*ga[i+1];
phase= pi/180.*(gp[i]*(1.-z)+z*gp[i+1]);
factor= complex(cos(phase),sin(phase));
g= mag*factor;
printf(" g magnitude, phase(rad) %le %le\n",mag,phase);
if(t>3.)// asymptotic results
	{
	sum=0.;
	for(i=0;i<10;i++) sum += (-.5*(sqrt(3.)-iu)*t*zai[i]).cexp()/aiz[i];
	f= sum*(-iu*pi/3.).cexp();
	for(i=0;i<10;i++) sum += (-.5*(sqrt(3.)-iu)*t*zaip[i]).cexp()/(zaip[i]*aipz[i]);
	g= sum;
	f.print(" f asymptotic=","");g.print(" g=","\n");
	mag= f.abs();phase= atan2( f.imaginary(),f.real());
	printf(" f mag=%le phase(rad)=%le\n",mag,phase);
	mag= g.abs();phase= atan2( g.imaginary(),g.real());
	printf(" g mag=%le phase(rad)=%le\n",mag,phase);
	return 2;
	}
else
	{
	mag= fa[i]*(1.-z)+z*fa[i+1];
	phase= pi/180.*(fp[i]*(1.-z)+z*fp[i+1]);
	factor= complex(cos(phase),sin(phase));
printf(" f magnitude, phase(rad) %le %le\n",mag,phase);
	f= mag*factor;
	}
return 0;
}


static double pr[51]={.04,.221,.377,.503,.602,.674,.723,
.754,.768,.771,.766,.756,.742,.726,.712,.699,.689,.682,
.679, .682,.69,.707,.732,.77,.827,.911,1.043,1.27,1.732,3.135,0.,
-2.522,-1.119,-.657,-.429,-.297,-.211,-.152,-.11,-.08,-.058,-.041,
-.029,-.019,-.012,-.008,-.004,-.001,.001,.002,.003};
static double pip[51]={.879,.84,.769,.678,.577,.469,.354,.265,.173,.091,.019,
-.043,-.113,-.139,-.174,-.202,-.224,-.24,-.251,-.257,-.260,-.26,-.256,-.252,
-.244,-.236,-.225,-.214,-.202,-.19,-.177,-.164,-.151,-.138,-.125,-.113,
-.101,-.09,-.08,-.07,-.061,-.053,-.045,-.039,-.032,-.027,-.023,-.018,-.014,
-.011,-.01};
static double qr[61]={-.135,-.314,-.458,-.568,-.646,-.694,-.717,-.718,
-.703,-.676,-.639,-.596,-.548,-.499,-.449,-.399,-.35,-.3,-.251,-.201,
-.15,-.096,-.035,.034,.119,.229,.385,.634,1.118,2.542,0.,-3.074,-1.65,
-1.166,-.918,-.762,-.654,-.573,-.508,-.454,-.409,-.369,-.333,-.301,-.273,
-.246,-.222,-.2,-.18,-.161,-.144,-.128,-.113,-.1,-.088,-.077,-.067,
-.058,-.05,-.042,-.033};
static double qi[61]={-.838,-.771,-.676,-.562,-.44,-.317,-.199,-.09,-.008
,-.094,-.166,-.226,-.274,-.311,-.338,-.357,-.368,-.372,-.371,-.365,
-.356,-.342,-.327,-.309,-.289,-.268,-.246,-.223,-.2,-.177,-.154,-.131,-.109,
-.088,-.067,-.048,-.031,-.014,-.0013,-.015,-.027,-.038,-.048,-.056,-.062,
-.068,-.072,-.075,-.078,-.079,-.079,-.079,-.078,-.077,-.075,-.072,-.070,
-.067,-.064,-.061,-.059};



// called p,q  or phat,qhat. various authors

int fockpq(double t, complex& p, complex& q)
{int i;double z,eps; complex factor,iu(0.,1.),sump,sumq;

if( t<= -3.)
	{
	z=t*t*t;factor=(-iu*(pi*.25+z/12.)).cexp();
	p= factor*(1.-iu*2./z+20./(z*z));
	q=-factor*(1.+iu*2./z-28./(z*z));
	return 1;
	}
if(  t >= 2.)
	{
	sump=0.; sumq=0.;
	for(i=0;i<10;i++) 
		{
		sump += (-.5*(sqrt(3.)-iu)*t*zai[i]).cexp()
		/(aiz[i]*aiz[i]);
		sumq +=  (-.5*(sqrt(3.)-iu)*t*zaip[i]).cexp()
		/(zaip[i]*aipz[i]*aipz[i]);
		}
	factor=-(-iu*pi/6.).cexp()/(2.*sqrt(pi));
	p= sump*factor;
	q= sump*factor;
	if(t>3.)return 1;
	}
eps=.00001;
i= (t+3.+eps)/.1; z= (t-i*.1)+3.;
//{printf("  z=%le t=%le i=%d\n",z,t,i);}
if(z<0.)z=0.;

if(t> -.1 && t< .1)
	{
	p=complex( .3-.5/(sqrt(pi)*t),pip[i]*(1.-z)+z*pip[i+1]);
	q=complex( -.25-.5/(sqrt(pi)*t),qi[i]*(1.-z)+z*qi[i+1]);
	}
if(t<=2.)p= complex(pr[i]*(1.-z)+z*pr[i+1],pip[i]*(1.-z)+z*pip[i+1]);
q= complex(qr[i]*(1.-z)+z*qr[i+1],qi[i]*(1.-z)+z*qi[i+1]);
return 0;
}




// lower case p*,q* scattering functions

int fockscat(double t, complex& p, complex& q)
{
double spi,spii,sqta,d;complex term,sum1,sum2,ee;int i;
spi=1./sqrt(pi);spii=spi*.5/t;
if(t>2.)//deep shadow
	{
	ee=(complex(0.,-5.*pi/6.)).cexp();
	for(i=0;i<5;i++)
		{
		d=aiz[i];
		term= (ee*t*zai[i]).cexp()/(d*d);
		sum1+=term;
		d=aipz[i];
		term=(ee*t*zaip[i]).cexp()/(d*d*zaip[i]);
		sum2+=term;
		}
	term= spi*(complex(0.,pi/6.)).cexp();
	p=spii-sum1*term;
	q=spii-sum2*term;
	}
else if(t<-3.) //brightly lit
	{ 
	sqta= sqrt( ABS(t) );
	term=((complex(0.,pi*.25)+complex(0.,1./12.)*t*t*t).cexp())*.5*sqta;
	p= spii +term*(1.+complex(0.,2.)/(t*t*t));
	q= spii -term*(1.-complex(0.,2.)/(t*t*t));
	}
else// interpolate -3<x<2
	{
	 fockpq(t,p,q);
	 d = .5/(sqrt(pi)*t);
	 q+= d;
	 if(abs(t)>=.1) p+=d;
	 else { spi= .314-.015*(t+.1)/.2;
			p=complex(spi,p.imaginary());
			spi= -.279+.026*(t+.1)/.2;
			q=complex(spi,q.imaginary());
			}/* interpolate real part p */
	}
return 0;
}

//radiation functions
// very crude- do not use if possible to avoid!!!!!!
int fockrad(double t, complex& G, complex& Gtwiddle)
{
if(t>3.)
	{G=1.8325*((complex(-.8823*t,-.5094*t)).cexp());//best t>4
	 Gtwiddle=complex(0.,0.);
	}
else if(t<-2.)
	{G=complex(2.,-.5/(t*t*t));// good for t<-2.5
	Gtwiddle=complex(2.,.5/(t*t*t)) *complex(0,-t);//best t<-4.
	}
else//interpolate- 
	//G plot only from -2 to 3 
	// f,g plots -3 to 4 still need -4 to -3 for Gtwiddle use above for G
	{
	G= complex( 2./( 1.+ .43* exp(t*1.65)), -t/(1.+2.*t*t));
	Gtwiddle=complex( .5/(1.+t*t) , 4.*exp( -.87*(t+2.)));

	}
return 0;
}


//coupling functions 
// u=  x^1.5  exp 3ipi/4 integral (w2'/w2)exp-jxt dt
// v= .5 x^.5 exp ipi/4  integral (w2/w2')exp-jxt dt 
// v1= x^1.5 exp 3ipi/4 integral t(w2/w2')exp-jxt dt  

fockc(double t, complex& u, complex& v, complex& v1)
{double sp,st,y;complex i,q,sum1,sum2,sum3,p;int k;
if(t<0.){q=complex(0.,0.);u=q;v=q;return 0;}
i=complex(0.,1.);
q=(complex(0.,pi*.25)).cexp();
p=(complex(0.,-pi/3.)).cexp();
sp=sqrt(pi);
st=sqrt(t);
y=t*st;
if(t<.6)
	{
	u= 1.+y*(-sp*.5*q+y*(5./12.*i+y*(sp*5./64./q-3.701e-3*y)));
	v= 1.+y*(-sp*.25*q+y*(7./60.*i+y*(sp*7./512./q-4.141e-3*y)));
	v1= 1.+y*(sp*.5*q+y*(-7./12.*i+y*(-sp*7./64.+4.555e-3*y)));
	}
else         //t>=.6
	{
	for(k=0;k<10;k++)
		{
//sum1.print(" sum1=","");
//sum2.print(" sum2=","\n");
		sum1+= (-i*t*p*zai[k]).cexp();
		sum3+= (-i*t*p*zaip[k]).cexp();
		sum2+= (-i*t*p*zaip[k]).cexp()/(p*zaip[k]);
		}
	u= q*sp*2.*sum1*y;
	v=sp*sum2/q*st;
	v1= q*sp*2.*sum3*y;
	}
return 0;
}


main()
{int i;double d;complex p,q,g,gt,u,v;
for (i=0;i<51;i++)
	{d= -3.+(i)*(6./50.);
	fockpq(d,u,v);
	printf(" x=%le\n",d);
	u.print(" phat=","");
	v.print(" qhat=","\n");
	}
u= complex(0.,-pi/4.); u=u.cexp();
for (i=0;i<51;i++)
	{d= -6.+(i)*(12./50.);
	fockscat(d,p,q);
	printf(" x=%le\n",d);
	p.print(" p= ","");
	q.print(" q=","\n");
	p*=u;q*=(u);
	p.print(" e*p= ","");
  (q).print(" e*q=","\n");
	}
for (i=0;i<21;i++)
	{d= -4.+(i)*(8./20.);
	printf(" x=%le ",d);
	fockfg(d,p,q);
	p.print(" f=","");
	q.print(" g=","\n");
	}


for (i=0;i<=20;i++)
	{d= +(i)*(.6/20.);
	printf(" x=%le\n",d);
	fockc(d,u,v,g);
	u.print("  u=","");
	v.print("  v=","\n");
	g.print("  v1=","\n");
	}
// x=.6 and x=.65 both will be sum not series
	{d= .65;
	printf(" x=%le ",d);
	fockc(d,u,v,g);
	u.print("  u=","");
	v.print("  v=","\n");
	g.print("  v1=","\n");
	}

for (i=0;i<20;i++)
	{d=.5 +(i)*((5.-.5)/20.);
	printf(" x=%le\n",d);
	fockc(d,u,v,g);
	u.print("  u=","");
	v.print("  v=","\n");
	g.print("  v1=","\n");
	}
/*
for (i=0;i<20;i++)
	{d= -6.+(i)*(12./20.);
	printf(" x=%le ",d);
	fockrad(d,g,gt);
	g.print(" G=","");
	gt.print(" Gt=","\n");
	}
*/
}

