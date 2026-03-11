/*  Generalized Hypergeometric Function and Meijer G
	Pade summation used for former
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
#include "Cmatrix.hpp"

#define maxterm 100
#define abstol 1.e-8
#define reltol 1.e-5
#define errorcode -1.e60
#define pi 3.141592653589793238462643383279
#define lim 15.
#define max(a,b) ((a)>(b)?(a):(b))
complex cgamma(complex&);

// generalized hypergeometric function pFq or pFq-1 
// if 0<=k<=q then the k-th factor (c)k in denominator omitted pFq-1
// if k=-1,say, then pFq= sum x^n/n! (a[0])n..../{(c[0])n.....}
// pade approx. version
int szused,sizep=60;

complex ghgF(int p,int q,complex x, Cvector& a, Cvector& c, int k)
{complex sum(1.,0.),oldsum,power(1.,0.),btm;double test,fact;
complex f,t,pc,qc,*h,hold;int j,kk; double diff,tk,z;
if( sizep%2)sizep++;//size must be even integer
h=new complex[sizep+1] ;

Cvector denom(max(q,1)),numer(max(p,1));
// use 2F1 if p=2,q=1
//x.print(" x in ghgf=","\n");
//printf(" in ghgF: p q %d %d k=%d\n",p,q,k);
fact=1.;
for( j=0;j<p;j++)
	{
	//(a.element(j)).print(" a=","\n");
	numer.setelement(j,a.element(j));
	//(numer.element(j)).print(" numer=","\n");
	}
for( j=0;j<q;j++) {denom.setelement(j,c.element(j));
		//(denom.element(j)).print(" denom=","\n");
		//c.element(j)).print(" a=","\n");
			}
//power.print(" power=","\n");

h[0]=1.;kk=0;hold=1.;
while(kk<sizep){
	kk++;
//test t[n]=x^n geometric series 1+x+x^2+...=1/(1-x)
//	t*=xx;
	btm=1.;
	for(j=0;j<p;j++){power *= numer.element(j);
			//power.print(" power=","\n");
			numer.setelement(j,numer.element(j)+1.);}
	for(j=0;j<q;j++){if(j!=k)btm *=denom.element(j);
			denom.setelement(j,denom.element(j)+1.);}
//btm.print(" btm=","\n");
power *= (x/(btm*fact));
//power.print(" power=","\n");
	t=power;
	h[kk]=h[kk-1]+t;
	fact+=1.;
	qc=0.;j=kk;
	do	{
		pc=h[j-1];
		if(h[j]==pc)
			{
			h[j-1]=-errorcode;
			if(pc==-errorcode)h[j-1]=qc;
			}
		else h[j-1]= qc+1./(h[j]-pc);
		qc=pc;j--;
		}while(j);
	if(!(kk%2))
		{
		diff=(hold-h[0]).abs();
		if(diff<abstol || diff<reltol * h[0].abs())break;
		hold=h[0];
		}
	}
szused=kk;
hold=h[0];
delete h;
return  hold;
}



complex Gmeijer(int m, int n, int p, int q,Cvector& a,Cvector& c, complex& x)
{/* Meijer G function*/
complex coef,sum,denom,differ;double re;int j;
if(m<1 ||m>q){fprintf(stderr," Meijer:illegal m\n");return errorcode;}
if(n<0 ||n>p){fprintf(stderr," Meijer:illegal n\n");return errorcode;}
Cvector aa(max(1,q)),cc(max(1,p));
if(p>q)    /* p>q implies for pFq-1 as p'Fq' p'>q'+1*/
	{
	//printf(" Meijer G: p>q\n");
	for(j=0;j<q;j++) aa.setelement(j, 1.-c.element(j));
	for(j=0;j<p;j++) cc.setelement(j, 1.-a.element(j));
	sum=Gmeijer(n,m,q,p,aa,cc,1./x);
	return sum;
	}
else if(p==q)
	{if( x.abs()>=1.){//printf(" p=q Meijer\n");
		return complex(errorcode,0.);}
	}
for(int k=0;k<m;k++)//arrays based at zero(default)
	{
	coef=complex(1.,0.);
	denom=complex(1.,0.);

for(int j=0;j<p;j++)
	{cc.setelement(j,1.+c.element(k)-a.element(j));
	//(cc.element(j)).print(" numer=","\n");
	}
for( j=0;j<q;j++) {aa.setelement(j,1.+c.element(k)-c.element(j));
		//(aa.element(j)).print(" denom=","\n");
			}

//coef.print(" coef=","\n");
	for( j=0;j<n;j++) if(j!=k)coef *= cgamma(1.-aa.element(j));
//coef.print(" coef=","\n");
	for( j=0;j<n;j++) coef *= cgamma(cc.element(j));
//coef.print(" coef=","\n");
	for( j=n;j<p;j++)denom*= cgamma(1.-cc.element(j));
//denom.print(" denom=","\n");
	for( j=m;j<q;j++) denom *= cgamma(aa.element(j));
	coef *= (x^c.element(k))/denom;
//denom.print(" denom=","\n");
//printf(" before call ghgF ");coef.print(" coef=","\n");
//sum.print(" sum=","\n");
	sum += coef*ghgF(p,q,((p-m-n)%2?-1.:1.)*x,cc,aa,k);
//printf(" after call ");sum.print(" sum=","\n");
	}
return sum;
}

main()
{
complex a,b,c;double x,y;
/*
while(1)
	{printf(" enter complex x\n");scanf("%le%le",&x,&y);
	if(x==0. && y==0.)break;
	a=complex(x,y); (cgamma(a)).print(" gamma=","\n");
	}
*/
Cvector A(10),B(10);
while(1)
	{printf(" enter complex x\n");scanf("%le%le",&x,&y);
	if(x==0. && y==0.)break;a=complex(x,y);
	A.setelement(0,complex(3.,0.) );
	b=Gmeijer(1,0,0,1,B,A,a); b.print(" answer="," ");
	((a^3.) * ((-a).cexp())).print("","\n");
	}

while(1)
	{printf(" 1F1:enter complex x\n");scanf("%le%le",&x,&y);
	if(x==0. && y==0.)break;b=complex(x,y);
	printf(" enter complex a\n");scanf("%le%le",&x,&y);a=complex( x,y);
	printf(" enter complex c\n");scanf("%le%le",&x,&y);c=complex( x,y);
	A.setelement(0,(1.-a));
	B.setelement(0,complex(0.,0.));
	B.setelement(1,(1.-c));
	b=Gmeijer(1,1,1,2,A,B,(-b))*cgamma(c)/cgamma(a);
	b.print(" 1F1(a,c,x) =","\n");
	}
}

error(char* msg,int index,int s)
	{
	printf("%s%d %d\n",msg,index,s);
	exit(1);
	return 0;/* keep compiler happy*/
	}

void Cvector::check( int index)
	{int loc;
	loc=index-(*this).base;
	if (loc<0)
		error(" Cvector error: index-base<0",index,(*this).base);
	if (loc>= (*this).size)
		error(" Cvector error: index too large",index,(*this).size);
	}

void Cmatrix::check( int i,int j)
	{int loc;
	loc=i-(*this).base;
	if (loc<0)
		error(" Cmatrix error: row index-base<0",i,(*this).base);
	if (loc>= (*this).rowkt)
		error(" Cmatrix error: row index too large",i,(*this).rowkt);

	loc=j-(*this).base;
	if (loc<0)
		error(" Cmatrix error: column index-base<0",j,(*this).base);
	if (loc>= (*this).colkt)
		error(" Cmatrix error: column index too large",j,(*this).colkt);

	}


void Cmatrix::init(int row, int col, int b)
	{
	if(row<1||col<1)error(" matrix needs at least one row/col ",row,col);
	rowkt=row;colkt=col;base=b;
	m= new complex *[rowkt];
	for(int j=0;j<rowkt;j++)m[j]= new complex[colkt];
	complex zero;complex one(1.,0.);
	//initialize to identity matrix. 
	for(int i=0;i<rowkt;i++){for(j=0;j<colkt;j++)m[i][j]=zero;
		m[i][i]= one;// omit this line for init. to zero matrix
		}
	}

Cmatrix::Cmatrix(Cmatrix& a)
	{ init(a.rowkt,a.colkt,a.base);
	 for(int i=0;i<a.rowkt;i++)for(int j=0;j<a.colkt;j++)m[i][j]=a.m[i][j];
	}

Cmatrix::~Cmatrix()
	{for(int i=0;i<rowkt;i++)delete m[i];
	delete m;
	}


complex& Cvector::element(int i)
	{
	check(i);return head[i-base];
	}

void Cvector::setelement( int i, complex& value)
	{
	check(i);
	head[i]=value;
	}
complex& Cmatrix::element(int i,int j)
	{
	check(i,j);return m[i-base][j-base];
	}

void Cmatrix::setelement( int i, int j, complex value)
	{
	check(i,j);m[i][j]=value;
	}


void Cmatrix::operator=(Cmatrix& a)
	{
	if(this==&a)return;
	for(int i=0;i<rowkt;i++)delete m[i];
	delete m;
	init(a.rowkt,a.colkt,a.base);
	for(i=0;i<a.rowkt;i++)for(int j=0;j<a.colkt;j++)m[i][j]=a.m[i][j];
	}

 Cvector::Cvector( Cvector& orig) //copy
	{
	size=orig.size;base=orig.base;head=new complex[orig.size];
	for(int i=0;i<size;i++)head[i]=orig.head[i];
	}

Cmatrix operator*(Cmatrix&a, Cmatrix& b)
	{
	int i,j,k;complex zero;
	if(a.colkt!=b.rowkt)
		error(" cannot multiply matrices colkt,rowkt ",a.colkt,b.rowkt);
	Cmatrix prod(a.rowkt,b.colkt,a.base);
	for(i=0;i<a.rowkt;i++)for(j=0;j<b.colkt;j++)
		{prod.m[i][j]=zero;
		for(k=0;k<a.colkt;k++)prod.m[i][j] += a.m[i][k]*b.m[k][j];
		}
	return prod;
	}	


Cvector operator*(Cmatrix& a, Cvector& b)
	{
	int k;complex zero;
	if(a.colkt!=b.size)error(" cannot mult. matrix,vector ",a.colkt,b.size);
	Cvector prod(b.size,b.base,zero);
	for(int i=0;i<a.rowkt;i++)
		{prod.head[i]=zero;
		for(k=0;k<a.colkt;k++)prod.head[i] += (a.m[i][k])*(b.head[k]);
		}
	return prod;
	}	

complex Cvector::operator*(Cvector& rvalue) //dot product
	{int i; complex *element,*relement; complex sum(0.,0.);// sum 
	if(rvalue.size!= size)
		error(" dot product unequal length vectors ",size,rvalue.size);
	for(i=0,element=head,relement=rvalue.head;i< size;
		i++,element++,relement++)
		sum += *element * *relement;
	return sum;
	}

void Cvector::operator=(  Cvector& rhs)
	{
	if( this== &rhs)return;
	if(size!=rhs.size)
		{
		delete head;
		head= new complex[rhs.size];
		}
	for(int j=0;j<size;j++){head[j]=rhs.head[j];}
	}
