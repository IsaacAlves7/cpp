/* prototypes for C Mathematical Function Handbook
from C Mathematical Function Handbook by Louis Baker
Copyright 1991 by Louis Baker. All rights reserved.
*/
#include "complex.h"

/* Chapter 3: powers and roots file:pr.c*/
double power( double x, double n);
double root(double x,int n);
double square_rt(double x);
double cube_rt(double x);
double Euclidd(double xin, double yin);

/* Chapter 3: roots of polynomial equations file:polyrt.c*/
solvq(struct complex *b,struct complex *c,struct complex *ans1,
	struct complex *ans2);
int cubic(double a,double b,double c,double *r1,double *r2,double *r3);
int quartic(double a,double b,double c,double d,
struct complex *r1,struct complex *r2,struct complex *r3,struct complex *r4);
ctreal(struct complex *x,double y,struct complex *ans);
cpow(struct complex *x,struct complex *y,struct complex *ans);
ccubic( struct complex *a1,struct complex *a2,struct complex *a3,
struct complex *r1,struct complex *r2,struct complex *r3);
cquartic( struct complex *a1,struct complex *a2,struct complex *a3,
struct complex *a4,
struct complex *r1,struct complex *r2,struct complex *r3,struct complex *r4);


/* Chapter 4: elementary functions file:elem.c*/
double ln(double x);
double expon(double x);
double sine(double xi);
double arc_tan(double xi);
double arc_sine(double xi);
double hyper_sin(double x);
double hyper_cos(double x);
double hyper_tan(double x);
double arc_hyper_sin(double x);
double arc_hyper_cos(double x);
double arc_hyper_tan(double x);
double cosine(double x);
double arc_cosine(double x);
double tangent(double x);
double tangnt(double x);
double arc_tangent(double y, double x);

/* Chapter 4: functions for complex analysis file: complex.c*/
double sign(double x);
double argmt(double  x,double y);
double Argmt(double x,double y);
clog( struct complex *x,struct complex *ans);
csqrt(struct complex *z,struct complex *ans);
polarxy(double r, double angle,	double *x, double *y);
cexp( struct complex *x,struct complex *ans);
ctrig( struct complex *z,struct complex *ccos,struct complex *csin);
csin(struct complex *x,struct complex *ans);
ccos(struct complex *x,struct complex *ans);
ctan(struct complex *x,struct complex *ans);
ccot(struct complex *x,struct complex *ans);
csinh(struct complex *x,struct complex *ans);
ccosh(struct complex *x,struct complex *ans);
casin(struct complex *x,struct complex *ans);
cacos(struct complex *x,struct complex *ans);
catan(struct complex *x,struct complex *ans);
printc(struct complex *z);

/* Chapter 5: exponential integral and relatives file: expi.c*/
double e1(double x);
double si(double x);
double ci(double x);
/*double en(double x);  en,f,g,eicf,fcf,f01,g01,g2,f2,ceii
not intended for direct user invocation
not made static for special cases where direct call desired */
double e1s(double x);
double en(double x, double n);
double eicf(double x);
double e(double x,int n);
double E1(double x);
double Eict(double x);
double Eis(double x);
double eias(double x);
double ei(double x);
double Ei(double x);
double li(double x);
double alpha(double x,int n);
double beta(double x,int n);
double c01(double x);
double s01(double x);
cei( double x,double y, double k,double toler,double *u,double *v,int *n);
ceii( double x,double y, double k,double toler,double *u,double *v,int *n);
cexpint(struct complex *z,double k,double toler,struct complex *ans,int *iter);
double f( double x);
double g( double x);
double f2( double x, double y);
double g2( double x, double y);
double f01( double x);
double g01( double x);
double fcf( double x,int n,  double a[], double b[]);
double c2(double x, double y);
double s2(double x, double y);

/* Chapter 6: Gamma function and relatives file: g.c*/
double gamma(double x);
double loggam(double x);
double P(double a,double x);
double incgam(double a, double x);
double BigGamma(double a, double x);
double Pgamma(double a, double x);
double SmallGamma(double a, double x);
double fac(double x);
double pochhammer(double z, int n);
double gammaqd(double xin);
cgamma(struct complex *zz,struct complex *ans, struct complex *loggam);
cdigamma(struct complex *x,struct complex *ans);
double polygam(double z,int n);

/* Chapter 6: Digamma function and relatives file: digam.c */
double digamma(double x);
double digam(double x);
double pg1(double x);
double pg2(double x);

/* Chapter 7: Error and relate functions. file: pdfs.c*/
cerror( struct complex *a, struct complex *b, double eps);
double erf( double z, double *erfc);
double dawson(double x);
fresnel(double z, double *fci, double *fsi);
pdisp( struct complex *zetai,double eps,
 struct complex *zeeo,struct complex *zeeprimo,int iter);
cfsmall(double eps);/* NOT intended for direct use by users*/
cfbig(double eps);/* NOT intended for direct use by users*/
double ritchie( double x);
/* Chapter 7: Complementary error function for complex arg. file: cerfc.c*/
cerfc( struct complex *z, struct complex *ans);
/* Chapter 7: Error and relate functions. file: ierfc.c*/
double ierfc(double z, int n);
double ierfcf(double x, int n);
ierfctable (double z, int n, double table[]);
/* Chapter 7: Error and relate functions. file: boehmer.c*/
double boehmer( double x, double nu, int type);
double Si(double x);
double Ci(double x);
double Shi(double x);
double Chi(double x);
double Fresnel( double x, int type);
/* ba: asympt. boehmer. not for direct use*/
double ba(double x,double nu, int type);

/* Chapter 8: Legendre Functions file: leg.c*/
double plm( int l,int m, double x);
double pl0(int l,double x);
double ql0(int l, double x);
double qlm(int l,int m,double x);
double pli(int l,double z) ;
double qli(int l,double z);
double legendrea(int m,int n,double x,int real);
qleg(int m,int nmax,double x,int real,double r[],double q[]);
double qnu(double x,double nu,int real);
double pnu(double x,double nu);
/* Chapter 8: Legendre Functions file: Mehler.c*/
double Mehler(double x,double z,int mm) ;
double Mehler0(double x,double z);
/* Chapter 8: Legendre Functions file: gaut.c*/
leg1(double x,int a,int nmax,double p[]);
leg2(double x,int m,int nmax,int d,double q[]) ;
leg3(double  x,int n,int mmax,int d,double q[]);
legend1(double x,double alpha,int nmax,int d,double p[]);
legend2(double x,double a,int m,int nmax,int d,double p[]);
conical(double x,double tau,int nmax,int d,double p[]);
toroidal(double x,int m,int nmax,int d,double q[]);
double arccosh(double x);
double conicalt(double theta,double tau);
double conicala(double x,double tau);
/* Chapter 8: Legendre Functions file: torp.c*/
double Ptoroidal(int n, double x);

/* Chapter 9: Bessel Functions file: bessr.c*/
double j0(double x);
double j1(double x);
double y0(double x);
double y1(double x);
double i0(double x);
double i1(double x);
double k0(double x);
double k1(double x);
void ke(double x, double *,double *,double *, double *);
void be(double x, double *, double *, double *, double *);
CTHET(double X, double *PARTR, double *PARTI);
CPHI(double X, double *PARTR, double *PARTI);
/* Chapter 9: Bessel Functions file: besst.c*/
double jn(double x, int n);
double yn(double x, int n);
double in(double x, int n);
double kn(double x, int n);
/* Chapter 9: Bessel Functions file: cbess.c*/
Bessel( int nn,struct complex *z,struct complex *j,struct complex *y,
	struct complex *h2,struct complex *jprime,
	struct complex *yprime,struct complex *h2prime);
bessel( int nn,struct complex *z,struct complex *j,struct complex *y,
	struct complex *h2,struct complex *jprime,
	struct complex *yprime,struct complex *h2prime,int *ivalck);
double jint(int n, double z);
ibess(struct complex *z,struct complex *i, int n);
kbess(struct complex *z,struct complex *k,int n);
kelvin(int n,double x,struct complex *be,struct complex *ke);
forward(struct complex *z,struct complex *ratio
	, int idim, struct complex *r1);
backward(struct complex *z,struct complex *ratio
	, int idim);
cbess(struct complex *z,struct complex *j0,struct complex *j1,
	struct complex *y0,struct complex *y1,struct complex *h20,
	struct complex *h21,int *ivchk);
/* Chapter 9: Bessel Function Zeros  file:bessz.c*/
double fi(double y);
double jass(double z,double nu);
besspq(double a, double x,double *pa,double *qa,double *pa1,double *qa1);
int zerobes(double a,int n, double z[],int d);
/* Chapter 10: spherical Bessel Functions  file:spbn.c*/
double sjn(double x, int n);
double syn(double x, int n);
double skn(double x, int n);
double sinb(double x, int n);
double si0(double x);
double si1(double x);
double sim1(double x);
double sim2(double x);
double sj0(double x);
double sj1(double x);
double sy0(double x);
double sy1(double x);
double sk0(double x);
double sk1(double x);

/* Chapter 10: Bessel Functions file:abb.c*/
/* bessel functions of order 1/3 and -1/3 */
double jt( double z);
double jmt( double z);
double it( double z);
double imt( double z);
/* bessel functions of real order*/
double jbes( double z, double nu);
double jas( double z, double nu);
double ybes( double z, double nu);
double kbes( double z, double nu);
double ibes( double z, double nu);
/* Airy functions, their derivatives and integrals*/
double smallf(double z);
double smallg(double z);
double ai(double z);
double bi(double z);
double dai(double z);
double dbi(double z);
double bigf(double z);
double bigg(double z);
double IAi(double z);
double IBi(double z);
/* Anger and Weber functions*/
aw( double nu, double z, double *jj, double *e);
/* Chapter 11: Integrals of Bessel Functions file:bessi.c*/
double ji(double x);
double jii(double x);
double jin(double x, int n);
double jiotn(double x, int n);
double ki(double x);
double kii(double x);
double ii0m1t(double x);
double ij0t(double x);
double iy0t(double x);
double ik0t(double x);
double ii0tas(double x);
double ik0tas(double x);
double y0i( double x);
double i0i( double x);
/* Chapter 11: Integrals of Bessel Functions file:bickley.c*/
double bii(double x);/* bii not for user direct call*/
double bickley(double x, double rr);
double Jrn(double x, double r, int n);
/* Chapter 11: Integrals of Bessel Functions file:simp.c*/
double adsimp(double a,double b,double eps,double (*f)());
double simp(double a,double da,
	double fa,double fm,double fb,double area,double est,
	double eps,double (*f)()); 
/* Chapter 12: Bessel Functions file:struv.c*/
double StruveH(double nu, double x);
double StruveL(double nu, double x);
/* Chapter 12: Bessel Functions file:struvl.c*/
double h0(double x);
double h0a(double x);
double h1(double x);
double h1a(double x);
double l0(double x);
double l0a(double x);
double l1(double x);
double l1a(double x);
/* Chapter 12: Integrals of Struve Bessel Functions file:struvi.c*/
double ModStruveI( double x);
double StruveI(double x);
double StruveIot(double x);
/* Chapter 12: Bessel Functions file:iaw.c*/
double gamtab(double x);
iaw(int m, struct complex *s, struct complex *ans);

/* Chapter 13: Confluent Hypergeometric Functions file: chf.c */
Cpow(struct complex *x,struct complex *y,struct complex *ans);
c1f1(struct complex  *a,struct complex  *c,struct complex  *x
	,int top,struct complex  *ans);
cu(struct complex  *a,struct complex  *c,struct complex  *x
	,struct complex  *ans);
/* Chapter 13: Confluent Hypergeometric Functions file: chfs.c */
Jbessel(struct complex *order,struct complex *arg,struct complex *ans);
Ibessel(struct complex *order,struct complex *arg,struct complex *ans);
Kbessel(struct complex *order,struct complex *arg,struct complex *ans);
Airy(struct complex *z, struct complex *ans);
BiAiry(struct complex *z, struct complex *ans);
bateman(struct complex *nu,struct complex *arg,struct complex *ans);
cunningham(struct complex *n,struct complex *m,
	struct complex *x,struct complex *ans);
toronto(struct complex *m,struct complex *n,
	struct complex *r,struct complex *ans);
charlier(int n,struct complex *nu,struct complex *x,struct complex *ans);
Laguerre(struct complex * a,int n,struct complex * x,struct complex *ans);
/* Chapter 13: Confluent Hypergeometric Functions file: chfw.c */
Mwhit(struct complex *k,struct complex *mu,
	struct complex *x,struct complex *ans);
Wwhit(struct complex *k,struct complex *mu,
	struct complex *x,struct complex *ans);
/* Chapter 13: Confluent Hypergeometric Functions file: u.c */
double uabx(double a, double b, double x, double eps, double *uprime);
/* following not intended for direct user call*/
brec(double a,double b,int k,double *f,double *g,double x);
chu(double a,double b,double x,int kmax,double eps,double u[],double *uprime);

/* Chapter 14:  Coulomb Wave Functions file: cwf.c */
/* only coulombf and coulombg are intended for direct user call*/
double coulombf(double eta,double rho,int l);
double coulombg(double eta,double rho,int l);
double gl(double q);
cwfa(int l, double eta,double rho, double *fl,double *gl, 
	double *flp, double *glp);

/* Chapter 15: Hypergeometric and relatives hyperg.c*/
double f21(double a,double b,double c,double x);
double f12(double a,double b,double c,double x);
double F01(double a,double x);
/* Chapter 15: Hypergeometric and relatives f211.c-f213.c*/
cpochhammer(struct complex *x, int n, struct complex *ans);
cf21( struct complex *a,struct complex *b,struct complex *c,
	struct complex *x,struct complex *ans);
f21big(struct complex *ain,struct complex *bin,struct complex *c,
	struct complex *x,struct complex *ans);
f211(struct complex *a,struct complex *b,struct complex *c,
	struct complex *x,struct complex *ans);

/* Chapter 15: Legendre P for complex parameters */
cp(struct complex *z,struct complex *mu,struct complex *nu,
      struct complex *ans);

/* Chapter 15: Legendre Q for complex parameters */
cq(struct complex *z,struct complex *mu,struct complex *nu,
      struct complex *ans);


/* Chapter 16:  Elliptic Functions file: cje.c */
cjef( struct complex *u, double m, struct complex *sn, 
	struct complex *cn, struct complex *dn);
double ratmp(double m);
double ratmm(double m);
double solvem(double ratio,double mu,double ml,double (*ratv)());
double getm(double g2,double g3,int d); 
weier(struct complex *z,double g2,double g3,struct complex *p,
	struct complex *pp,double *mp,double *kp,double *ep,
	double *omegap,double *eta,double *ehp);
ctld(struct complex *v,struct complex *q,struct complex *t1,
	struct complex *t2,struct complex *t3,struct complex *t4);
sigma(struct complex *z,struct complex *ans,double k,double omega,
	double eta,double m,double g2,double g3);
zetaw(struct complex *z,struct complex *ans, double k,double omega,
	double eta,double m,double g2,double g3);
/* Chapter 16:  Elliptic Functions file: ct5.c */
ctheta(struct complex *v,struct complex *q,struct complex *ct1,
	struct complex *ct2,struct complex *ct3,struct complex *ct4);
double q(double m);
double mq(double q);
emf(struct complex *q,struct complex *m) ;
emft(struct complex *t,struct complex *m) ;
amc(struct complex *x,double m,struct complex *ans) ;
/* Chapter 17:  Elliptic Integrals file: te.c */
tek(int id, double m, double *k, double *e);
tef(double phi, double m, double sig, double *f, double *e);
double jzeta(double phi, double m);
double am(double u, double m);
jef(double u, double m,double *sn,double *cn,double *dn);
double heuman(double phi,double m);
theta(double v,double m,double *t1,double *t2,double *t3,double *t4);
neville(double  u,double m,double *ts,double *tc,double *td,double *tn);
/* Chapter 17:  Elliptic Integrals file: E3.c */
double e3(double n,double phi,double m,double sig);
/* Chapter 17:  Elliptic Integrals file: cef.c */
cef(double m,struct complex *z,struct complex *e,struct complex *f,double sig);
czeta(double m,struct complex *z,struct complex *ans,double sig);
/* Chapter 17:  Elliptic Integrals file: ce.c */
ce(struct complex *u,double m,struct complex *ce);

/* Chapter 18:  Weierstrass Functions file: invw.c*/
invp(double g2,double g3,struct complex *z,struct complex *ans);
/* Chapter 18:  Elliptic Functions Brent root finder file: brent.c */
double brent(double a,double b,double eta,double t, double (*f)());

/* Chapter 19:  Parabolic Cylinder Functions file:  */
double wpcf(double a, double xx);
double wa(double a, double xx);
double wairy(double a, double xx);
double upcf(double nu, double x);
double vpcf(double nu, double x);
double mwhit(double k, double mu, double z);
double wwhit(double k, double mu, double z);
double m(double a, double b, double z);
double u(double a, double b, double z);

/* Chapter 20:  Mathieu Functions file: mathieu.c */
int matheign(double q, int r, int odd,double *eigenv);
double mathieu(double x,double q, int r, 
	double eigenv[],int sol,int fnc, int norm);
int tmofa(double alfa, double *tm, double *dtm);
int coef();
double fj(int n);
double fy(int n);
double dy(int n);
double dj(int n);
double ds(int n);
double dc(int n);
double dds(int k);
double ddc(int k);
double ps(int k);
double pc(int k);
double dps(int k);
double dpc(int k);
sum(double (*func)());
bessinit(int sol,int n);
int math(double xx,double qq,int r, double cv,
	int sol,int fnc,int norm,double f[],int k[]);
int bounds(int k,double approx,double tola,double *cv,int coln);
int mfitr8(double tola,double *cv,double *dcv);
int mfcval(int n,int r,double qq,double *cv,int coln);

/* Chapter 21:  Spheroidal Wave Functions file: spherr.c */
double simb(double x,int n);
double dssph(double c, double s);
double fs(double c, double s);
double es(double c, double s);
double gs(double c, double s);
int ste(double d[],double e[],int n);
int solve(double c,int prolate,int odd,int order,int nmax);
int figi(double above[],double below[],int nmax);
double nmn(int m,int n,int limit);
double rho(int m,int n, int limit);
double betsph(int m,int r);
double gamsph(int m,int r);
double asph(int m,int r);
double bsph(int m,int r);
double csph(int m,int r);
int tridi(double bl[],double diag[],double ab[],double c[],int n);
int setd(int mm,int n);
double angular(int mm,int n,int kind,double eta,int imag);
double radswf(int m,int n,int kind,double eta);
double radial(int m,int n,int kind,double xi);
double sphjoin(int m,int n,int kind,int limit);
/* Chapter 22: Legendre Functions file: legp.c*/
legptable(double x, int ntop, int m);

/* Chapter 22:  Orthogonal  file: orthop.c */
/* backp not intended for direct user call*/
double backp(int n,int cc,int typep,double x);
double Pjacobi(double alpha,double beta,int n,double x) ;
double laguerre(double alpha,int n,double x);
double Cgegenbauer(double alpha,int n,double x);
double Tcheby( int n, double x);
double Ucheby( int n, double x);
double Plegendre(int n, double x);
double Hermite(int n, double x);

/* Chapter 23: Zeta (real arg.) and relatives */
double zeta( double r);
double bernoulli (int n);
double euler (int n);
double bernpoly ( int n, double x);
double eulerpoly ( int n, double x);
double debye (double x, int n);
double dilog(double x);
double clausen(double theta);
double zeta1( double z);
double zeta2( double z, double a);
double bigphi(double z, double s, double v);
double betacat( double s);
double lambda(double s);
double eta( double s);
double ifermi(double mu, double s);

/* Chapter 24:Stirling Numbers. file: stirl.c */
double stirl1( int m, int n);
double stirl2( int m, int n);
double stirlingf( int m, int n);
/* Chapter 24: Fibonacci numbers. file:fib.c */
double fib( int n);
/* Chapter 24: binomial coefficient. file:binom.c */
double binom(double n, double m);
/* Chapter 25: Numerical Methods: no functions*/

/*Chapter 26: Statistics. file:stath.c*/
/* independent module. No prototypes*/

/*Chapter 26: Statistics. file:random.c*/
/* random variates */
double erlang( int k, double mean);
double expon(double mean);
double cauchy();
double logistic(double a, double k);
int randi(int nvalue, float probd[]);
double uniform(double a, double b);
double normal( double mean, double sd, double *s1, double *s2);
double norm(double mean, double sd);
/* random number generators*/
double u32(void);
double u16(void);
double randm(void);
/* alternative exponential, cauchy, normal */
double ex(void);
double ca(void);
double na(void);

/*Chapter 26: Statistics.  file:dist.c*/
/*distributions*/
double binomial_dist(int x,int n, double p);
double neg_binomial_dist(int x,int y, double p);
double poisson_dist(int x,double lambda);
double hypergeometric_dist(int N, int m, int k, int x);
double erlang_dist(double b,int c,double x);
/* random numbers*/
int binomial(int n, double p);
int neg_binomial(int x, double p);
double chisq(int v);
double extreme (double a, double b);
int geometric(double p);
int hypergeometric(int N, int X, int smalln);
double lognormal(double median, double sigma);/* sigma=shape parameter*/
double pareto( double c);
int poisson(double lambda);
double weibull(double b, double c);
double binomial_coef(int n, int m);

/*Chapter 27: Misc.
	see also ritchie in pdfs.c,debye in ebznew.c*/
/*Chapter27: sievert integral file:misc.c*/
double sievert( double x, double theta);
/*Chapter27: Clebsh-Gordon coeff. & relatives file:wigner.c*/
int notint(double x);
int nothint( double x);
int triangle( double j1, double j2, double j);
double m1e(double x);
double wigner(double j,double j1,double j2,double m,
	double m1,double m2);
double ClebshGordon(double j1,double m1,double j2,double m2,
	double j,double m);
double wigner3j(double j,double j1,double j2,double m,double m1,
	double m2);
double delta(double a,double b,double c);
double gammn(double x);
double racah(double a,double b,double c,double d,double e,double f);
double wigner6j(double a,double b,double c,double d,double e,double f);
double Wigner6j(double a,double b,double c,double d,double e,double f);
double Wigner9j(double a,double b,double c,double d,double e,double f
	,double g,double h,double i);
double X(double a,double b,double c,double d,double e,double f
	,double g,double h,double i);
double V(double a,double b,double c,double A,double B,double C);
double CG(double j1,double j2,double j,double m1,double m2,double m);
/*Chapter27: dilogarithm function (see also dilog in ebznew.c) 
						file: dilog.c*/
double Dlog(double x);
/*Chapter27: "Abramowitz" functions file: F123.c*/
double af1s(double x);
double af2s(double x);
double af3s(double x);
double afa(double x,double n);
double af1(double x);
double af2(double x);
double af3(double x);
double af(double x,int n);
/*Chapter27: polylogarithm function file: polylog.c*/
double polylog(double x, double n);
/*Chapter27: Lobashevsky function NOT IN A&S file:Lob.c*/
double Lob(double x);
/*Chapter27: mu, nu and relatives function NOT IN A&Sfile: nu.c*/
double mu(double x,double ain,double bin);
double nu(double x);
double Nu(double x,double a);
double nui(double t);
/*Chapter 28:Scales no functions*/

/*Chapter 29: C++ no prototypes here*/
