/* 
Bessel functions for complex arguments
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.

bessel returns bessel functions and derivatives for complex arguments,
	integer orders.
jint(int n, double z) Jn(z)
ibess	In(z) complex z integer n
kbess	Kn(z) complex z integer n
kelvin	kelvin functions of order n, real argument

service routines (not intended for direct user call):

forward	forward recursion in order
backward backward recursion in order
cbess   lowest order functions calculated

based upon FORTRAN subroutine BESSEL by R. C. Lindberg 
published in Math Note 1 by Air Force Weapons Lab.
Kirtland AFB, NM (now Phillips Lab.)

*/

#include <stdlib.h>
#include <stdio.h>
#include "cmlib.h"
#include "protom.h"

static unsigned int btm,top;

int Bessel( int nn,struct complex *z,struct complex *j,
	struct complex *y,
	struct complex *h2,struct complex *jprime,
	struct complex *yprime,struct complex *h2prime)
{int ivalck,rv;
rv=bessel(nn,z,j,y,h2,jprime,yprime,h2prime,&ivalck);
return ivalck;
}

forward(z,ratio,idim,r1) int idim; struct complex *z,*ratio,*r1;
{
struct complex p,twooz,two,c1;
int i,imax;
two.x=2.;two.y=0.;
c1.x=1.;c1.y=0.;
imax=idim-1;
CDIV(twooz, two,(*z));
CASSN((ratio[0]),r1);
for (i=0;i<imax;i++)
	{
/*printf(" i,imax=%d %d\n",i,imax);
if( &(ratio[i])>top)printf(" femory top %x %x %d\n",ratio[i],top,i);
if( &(ratio[i+1])>top)printf(" fmem+1 top %x %x %d\n",ratio[i],top,i);
if( &(ratio[i])<btm)printf(" femory btm %x %x %d\n",ratio[i],btm,i);
if( &(ratio[i+1])<btm)printf(" fem+1 btm %x %x %d\n",ratio[i],btm,i);
*/
	CMULT(p,(ratio[i]),twooz );
	CTREAL(p,p,(double)(i+1));
	CSUB(p,p,c1);
	CDIV((ratio[i+1]),p,(ratio[i]));
	}
return 0;
}

#define limit 1.e-8
#define big   ( - ( errorcode ) )

backward(z,ratio,idim) int idim; struct complex *z,*ratio;
{
struct complex p,denom,twooz,two,c1;
int j,i;
two.x=2.;two.y=0.;
c1.x=1.;c1.y=0.;
CDIV(twooz, two,(*z) );
for (j=0;j<idim;j++)
	{CMPLX(ratio[j],0.,0.);}
CTREAL(denom,*z,.5/(idim));
CLET((ratio[idim-1]),denom);
for(i=idim-2;i>=0;i--)
	{
/*
printf(" backward: i=%d \n",i);
if( &(ratio[i])>top)printf(" bemory top %x %x %d\n",ratio[i],top,i);
if( &(ratio[i+1])>top)printf(" bmem+1 top %x %x %d\n",ratio[i],top,i);
if( &(ratio[i])<btm)printf(" bemory btm %x %x %d\n",ratio[i],btm,i);
if( &(ratio[i+1])<btm)printf(" bmem+1 btm %x %x %d\n",ratio[i],btm,i);
*/
	CTREAL(p,twooz, (double)(i+1) );
	CSUB(denom ,p ,ratio[i+1]);
	if( cabs(denom)>limit)
		{
		CDIV((ratio[i]),c1,denom);
		}
	else
		{
		CMPLX((ratio[i]),big,0. );
		}
	}
return 0;
}

bessel(nn,z,j,y,h2,jprime,yprime,h2prime,ivalck)
struct complex *z,*j,*y,*h2,*jprime,*yprime,*h2prime; int nn,*ivalck;
{
struct complex p,q,ci,c1,c0,cb,jz,j1,yz,y1,hz,h1,r1,aux1,aux2;
struct complex *b,*fy,*fh,yzadd,jzadd,hzadd,noverz;
/*struct complex b[1000];*/
int i,iret,check,idim,n;
n=nn;
CMPLX(ci,0.,1.);
CMPLX(c0,0.,0.);
CMPLX(c1,1.,0.);
/*CMPLX(cb,1.e10,0.);*/
CMPLX(cb, abs(errorcode) ,0.);
iret=0; *ivalck=0;
if(n<0)
	{
	n=-n;
	iret=1;
	}
if(n<=1)
	{
	cbess(z,&jz,&j1,&yz,&y1,&hz,&h1,&check);
	if(check)*ivalck=1;
	if(n==1)
		{
		CSET(j,j1);
		CSET(y,y1);
		CSET(h2,h1);
		if(cabs(*z)!=0.)
			{
			CDIV(p,c1,*z);
			CMULT(q,p,j1);
			CSUB(*jprime,jz,q);
			CMULT(q,p,y1);
			CSUB(*yprime,yz,q);
			CMULT(q,p,h1);
			CSUB(*h2prime,hz,q);
			if(iret)
				{
				if( !(n-((n>>1)<<1)) )return *ivalck;
				CTREAL((*j),(*j),-1.);
				CTREAL((*y),(*y),-1.);
				CTREAL((*h2),(*h2),-1.);
				CTREAL((*jprime),(*jprime),-1.);
				CTREAL((*yprime),(*yprime),-1.);
				CTREAL((*h2prime),(*h2prime),-1.);
				}
			return *ivalck;
			}/* z=0*/
		CMPLX(*jprime,0.,0.);
		CMPLX(*yprime,0.,0.);
		CMPLX(*h2prime,0.,0.);
		return *ivalck;
		}
	/*n=0*/
	CSET(j,jz);
	CSET(y,yz);
	CSET(h2,hz);
	CSUB(*jprime,c0,j1);
	CSUB(*yprime,c0,y1);
	CSUB(*h2prime,c0,h1);
	return *ivalck;
	}
if(cabs(*z)==0.)
	{
	CSET(j,c0);
	CSET(jprime,c0);
	CSET(y,cb);CTREAL(*y,*y,-1.);
	CSET(h2,cb);
	CSET(yprime,cb);
	CSET(h2prime,cb);CTREAL(*h2prime,*h2prime,-1.)
	/* *ivalck=1; default to 0 */
	return *ivalck;
	}
/*else*/
CDIV(p,c1,*z);
CTREAL(noverz,p,((double)n));
cbess(z,&jz,&j1,&yz,&y1,&hz,&h1,&check);
CSET(j,jz);
/* set up via malloc b array of idim length*/
idim=20*n;
b= malloc( sizeof( struct complex) *idim+50*8);
if(!b)
	{
	fprintf(stderr," bessel: malloc memory error backward\n");
	return (*ivalck=2);
	}
/*
btm=b;
top=btm+idim;
*/
backward(z,b,idim);
for(i=0;i<n;i++)
	{
	CMULT(p,*j,(b[i]));
	CSET(j,p);
	}
CSUB(q,noverz,(b[n]));
CMULT((*jprime),q,(*j));
CDIV(r1,y1,yz);
idim=n+2;
fy= malloc( sizeof(struct complex) *idim+50*8);
/*
btm=fy;
top=btm+idim;
*/
if(!fy)
	{
	fprintf(stderr," bessel: malloc memory error forward\n");
	return (*ivalck=3);
	}
forward(z,fy,idim,&r1);
CLET((*y),yz);
for(i=0;i<n;i++)
	{
	CMULT(p,(*y),(fy[i]));
	CSET(y,p);
	}
CSUB(q,noverz,fy[n]);
CMULT(*yprime,q,*y);
CDIV(r1,y1,yz);
CMULT(jzadd,*j,b[n]);
CMULT(yzadd,*y,fy[n]);
if( (*z).y ==0.)
	{
	CMULT(q,ci,*y);
	CSUB(*h2,*j,q);
	CMULT(q,yz,fy[n]);
	CMULT(aux2,q,ci);
	CMULT(p,jz,b[n]);
	CSUB(q,p,aux2);
	CMULT(aux1,noverz,*h2);
	CSUB((*h2prime),aux1,q);
	CSUB(hzadd,p,aux2);
	free(b);free(fy);
	}
else
	{
	CDIV(r1,h1,hz);
	idim=n+2;
	fh= malloc( sizeof(struct complex) *idim+8*50);
	if(!fh)
		{
		fprintf(stderr," bessel: memory error forward2\n");
		return (*ivalck=4);
		}
/*
btm=fh;
top=btm+idim;
*/
	forward(z,fh,idim,&r1);
	CSET(h2,hz);
	for(i=0;i<n;i++)
		{
		CMULT(q,*h2,fh[i]);
		CSET(h2,q);
		}
	CMULT(p,noverz,(*h2));
	CMULT(q,(*h2),fh[n]);
	CSUB((*h2prime),p,q);
	CMULT(hzadd,(*h2),fh[n]);
	free(fh);
	if(iret)
		{
		/*printf(" processing for negative n\n");*/
		if( !(n-((n>>1)<<1)) )return 0;
		CTREAL((*j),(*j),-1.);
		CTREAL((*y),(*y),-1.);
		CTREAL((*h2),(*h2),-1.);
		CADD(q,noverz,b[n]);
		CMULT((*jprime),q,(*j));
		CADD(q,noverz,fy[n]);
		CMULT((*yprime),q,(*y));
		CMULT(q,ci,(*yprime));
		CSUB((*h2prime),(*j),q);		
		}
/*printf(" fini bessel\n");*/
	}
/* diff eq. chk*/
return *ivalck;
}

struct complex wronsk;

cbess(z,j0,j1,y0,y1,h20,h21,ivchk)
struct complex *z,*j0,*j1,*y0,*y1,*h20,*h21;int *ivchk;
{
struct complex p,q,zsq,fact,zfact,j0add,j1add,ez,aux1,aux2,aux,cosp,sinp,
cii,c1,ci;
double skinv, fk,fn,fkfact,fkfi,tol=1.e-8,u,tolz=1.e-8,tolw=1.e-8;
int phase;
/*printc(z);printf(" = z entered cbess\n");*/
*ivchk=0;
CMPLX(cii,1.,1.);
CMPLX(c1,1.,0.);
CMPLX(ci,0.,1.);
if( cabs(*z)<= 10.)
	{
	CMPLX(*j0,1.,0.);
	CMPLX(*j1,1.,0.);
	CMPLX(*y0,0.,0.);
	CMPLX(*y1,1.,0.);
	CMULT(p,*z,*z);
	CTREAL(zsq,p,-.25);
	fk=1.;
	fkfact=1.;
	CMPLX(zfact,1.,0.);
	skinv=1.;
	do
		{/*series*/
		CMULT(p,zfact,zsq);
		CLET(zfact,p);
		fkfi=1./fkfact;
		CTREAL(fact,zfact,fkfi);
		CTREAL(j0add,fact,fkfi);
		fk+=1.;
		fkfact*=fk;
		fkfi=1./fkfact;
		CTREAL(j1add,fact,fkfi);
		CADD(*j0,*j0,j0add);
		CADD(*j1,*j1,j1add);
		CTREAL(p,j0add,skinv);
		CADD(*y0,*y0,p);
		CTREAL(p,j1add, (skinv+skinv+1./fk));
		CADD(*y1,*y1,p);
		skinv+=1./fk;
		}while(cabs(j1add)>cabs(*j0)*tol || cabs(j1add)>cabs(*j1)*tol);
	CMULT(p,*j1,*z);
	CTREAL(*j1,p,.5);
	/*printc(j0);printc(j1);printf(" j0, j1 after series\n");*/
	if(cabs(*z)< tolz)
		{
		CMPLX(*y0,-big,0.)
		CMPLX(*y1,-big,0.);
		CMPLX(*h20,big,0.);
		CMPLX(*h21,big,0.);
		*ivchk=0;
		return 0;
		}

	CTREAL(p,*z,.5);
	clog(&p,&q);q.x+=.5772156649;
	fkfi=1./1.570796326795;
	CMULT(p,q,*j0);
	CSUB(*y0,p,*y0);
	CTREAL(*y0,*y0,fkfi);
	CMULT(p,q,*j1);
	CMULT(aux1,*y1,zsq);
	CSUB(aux2,c1,aux1);
	CDIV(aux1,aux2,*z);
	CSUB(q,p,aux1);
	CTREAL(*y1,q,fkfi);
	CMULT(aux1,ci,*y0);
	CSUB(*h20,*j0,aux1);
	CMULT(aux1,ci,*y1);
	CSUB(*h21,*j1,aux1);

	goto wronskian;
	}
/*else*/
if( abs(z->y)<1.e-9 && z->x < 0.)
	{
	phase=1;
	z->x = -(z->x);
	}
else
	phase=0;
/* above added to handle phase= pi, ie z real <0 */
CTREAL(ez,*z,8.);
CTREAL(fact,*z,pi);
csqrt(&fact,&aux);
CDIV(fact,c1,aux);
ctrig(z,&aux1,&aux2);
CMULT(cosp,aux1,fact);
CMULT(sinp,aux2,fact);
CMULT(aux1,*z,ci);
CTREAL(aux1,aux1,-1.);
cexp(&aux1,&aux);
CMULT(aux1,aux,fact);
CMULT(zfact,aux1,cii);
u=0.;
retry:
fn=1.;
fk=1.;
CMPLX(p,1.,0.);
CMPLX(aux, (u-1.),0.);
CDIV(q,aux,ez);
CLET(fact,q);
do
	{
	fn+=2.;
	fk+=1.;
	CTREAL(aux,fact, -(u-fn*fn)/fk);
	CDIV(fact,aux,ez);
	CADD(p,p,fact);
	fn+=2.;
	fk+=1.;
	CTREAL(aux,fact, (u-fn*fn)/fk);
	CDIV(fact,aux,ez);
	CADD(q,q,fact);
	}while(fk<21. && cabs(fact)>cabs(q)*tol);
if( u <=0.)
	{
	CADD(aux1,p,q);
	CSUB(aux2,p,q);
	CMULT(aux,aux1,cosp);
	CMULT(fact,aux2,sinp);
	CADD(*j0,aux,fact);
	CMULT(aux,aux1,sinp);
	CMULT(fact,aux2,cosp);
	CSUB(*y0,aux,fact);
	CMULT(aux,ci,q);
	CSUB(aux,p,aux);
	CMULT(*h20,aux,zfact);
	u=4.;
	goto retry;
	}
CADD(aux1,p,q);
CSUB(aux2,p,q);
CMULT(aux,aux1,sinp);
CMULT(fact,aux2,cosp);
CSUB(*j1,aux,fact);
CMULT(aux,aux1,cosp);
CMULT(fact,aux2,sinp);
CADD(*y1,aux,fact);
CTREAL(*y1,*y1,-1.);
CMULT(aux,ci,q);
CADD(aux,p,aux);
CMULT(*h21,aux,zfact);
if(phase)
	{
	z->x = -(z->x);
	CMULT(p,*j0,ci);
	CTREAL(p,p,2.);
	CADD( *y0,*y0,p);
	CMULT(p,*y0,ci);
	CSUB(*h20,*j0,p);
	CTREAL(*j1,*j1,-1.);
	CMULT(p,*j1,ci);
	CTREAL(p,p,2.);
	CSUB(*y1,p,*y1);
	CMULT(p,*y1,ci);
	CSUB(*h21,*j1,p);
	}
wronskian:
/*wronskian*/
if( cabs(*z)<tolz) return 0;
CDIV(aux,c1,*z);
CTREAL(aux,aux,2./3.14159265358979);
CMULT(p,*j1,*y0);
CMULT(q,*j0,*y1);
CADD(q,q,aux);
CDIV(wronsk,p,q);
CSUB(wronsk,wronsk,c1);/* should be zero*/
if(cabs(wronsk)>tolw) *ivchk=10;
return *ivchk;
}


ibess(z,i,n) struct complex *z,*i;int n;
{
struct complex ci,j,zz,dummy,mult;
double theta,argmt();   int ck;
CMPLX(ci,0.,1.);
theta=argmt(z->y,z->x);
CMULT(zz,(*z),ci);
bessel(n,&zz,&j,&dummy,&dummy,&dummy,&dummy,&dummy,&ck);
if(ck)fprintf(stderr," ibess j= %le %le ck=%d\n",j.x,j.y,ck);
if(ck || abs(j.x)==abs(errorcode))
	{
	i->x=i->y=errorcode;
	return 1;
	}
CSET(i,j);
if(theta <= pi*.5 ||theta> pi)
	{
	CTREAL(dummy,ci, -pi*.5*n);
	cexp(&dummy,&mult);
	CMULT(*i,*i,mult);
	return 0;
	}
/*else*/
CTREAL(dummy,ci, pi*1.5*n);
cexp(&dummy,&mult);
CMULT(*i,*i,mult);
return 0;

}

kbess(z,k,n) struct complex *z,*k;int n;
{
struct complex ci,j,y,h1,h2,zz,mult,dummy;
double theta,argmt();   int ck;
CMPLX(ci,0.,1.);
/*printf("entered kbess %d\n",n);*/
theta=argmt(z->y,z->x);
if(theta <= pi*.5 ||theta> pi)
	{
	CMULT(zz,(*z),ci);
	/*printf(" calling bessel\n");*/
	bessel(n,&zz,&j,&y,&h2,&dummy,&dummy,&dummy,&ck);
	/*printf(" leaving bessel y=%le %le\n",y.x,y.y);*/
	if(ck || abs(y.x)==abs(errorcode))
		{
		k->x=k->y=errorcode;
		return 1;
		}
	CMULT(dummy,ci,y);
	CADD(h1,j,dummy);
	CTREAL(dummy,ci, pi*.5*n);
	cexp(&dummy,&mult);
	CMULT(dummy,h1,mult);
	CMULT(*k,dummy,ci);
	CTREAL(*k,*k, pi*.5);
	/*printf("leaving kbess\n");*/
	return 0;
	}
/*else*/
CMULT(zz,(*z),ci);
CTREAL(zz,zz,-1.);
bessel(n,&zz,&j,&y,&h2,&dummy,&dummy,&dummy,&ck);
if( ck || abs(h2.x)==abs(errorcode))
	{
	k->x=k->y=errorcode;
	return 1;
	}
CMULT(*k,h2,ci);
CTREAL(dummy,ci, -pi*.5*n);
cexp(&dummy,&mult);
CMULT(dummy,h2,mult);
CMULT(*k,dummy,ci);
CTREAL(*k,*k,-pi*.5);
return 0;
}

double jint(n,z) double z;int n;
{
struct complex x,dummy,ans;
int ck;
x.x=z;x.y=0.;
bessel(n,&x,&ans,&dummy,&dummy,&dummy,&dummy,&dummy,&ck);
return ans.x;
}

kelvin(n,x,be,ke) double x; struct complex *be,*ke;int n;
{
struct complex mult,h2,dummy,ans,ci,y,coeff,h1;
int ck;
/*CMPLX(mult,0., .75*pi);
cexp(&mult,&coeff);
printf(" exp 3/4pii=");printc(&coeff);printf("\n");
*/
CMPLX(coeff,-.707106781,.707106781);
CMPLX(ci,0.,1.);
CTREAL(ans ,coeff,x);
bessel(n,&ans,be,&y,&h2,&dummy,&dummy,&dummy,&ck);
if( ck || abs(y.x)==abs(errorcode))
	{
	be->x=be->y=ke->x=ke->y=errorcode;
	return 1;
	}
CMULT(dummy,ci,y);
CADD(h1,(*be),dummy);
CMULT(dummy,ci,h1);
CTREAL((*ke),dummy,pi*.5);
return 0;
}
