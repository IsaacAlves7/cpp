/* COMPLEX.H header file
 * use for complex arithmetic in C
 * see MULLER.C for support functions such as
	csqrt(),clog(),cexp(),argmt(),polarxy()

Copyright 1991 by Louis Baker. All rights reserved.
Based upon COMPLEX.H header file in 
Handbook of C tools for Scientists and Engineers by L. Baker

*/

#ifndef complex

#include <math.h>

#else

double sqrt();
struct complex { double x;
		 double y;} ;
#endif

static double TP,T2,T3;/* dummy static so no conflict with other modules*/

/* for below, X,Y are complex structures, and one is returned*/

#define CMULTR(X,Y) ((X).x*(Y).x-(X).y*(Y).y)
#define CMULTI(X,Y)  ((X).y*(Y).x+(X).x *(Y).y)
/* CMLT uses 1 fewer multiplication but more adds than CMULT not generally
recommended */
#define CMLT(Z,X,Y) {TP=(X.x+X.y)*(Y.x+Y.y);T2=X.x*Y.x;T3=X.y*Y.y;Z.y=TP-T2-T3;z.X=T2-T3;}

#define CDRN(X,Y)  ((X).x*(Y).x+(Y).y*(X).y)
#define CDIN(X,Y)  ((X).y*(Y).x-(X).x*(Y).y)
#define CNORM(X) ((X).x*(X).x+(X).y*(X).y)
/*#define CNRM(X) (X->x*X->x+X->y*X->y)
*/
#define CDIV(z,nu,de) {TP=CNORM(de);(z).x=CDRN(nu,de)/TP;(z).y=CDIN(nu,de)/TP;}
#define CONJG(z,X) {(z).x=(X).x;(z).y=-(X).y;}
/*#define CONJ(X) {(X).y=-(X).y}
*/
#define CMULT(z,X,Y) {(z).x=CMULTR((X),(Y)); (z).y=CMULTI((X),(Y));}
#define CADD(z,X,Y) {(z).x=(X).x+(Y).x;(z).y=(X).y+(Y).y;}
#define CSUB(z,X,Y) {(z).x=(X).x-(Y).x;(z).y=(X).y-(Y).y;}
#define CLET(to,from) {(to).x=(from).x;(to).y=(from).y;}
/* For a "bulletproof" form of cabs use 
	#define cabs(X) Euclidd((X).x,(X).y)
*/
#ifndef cabs
#define cabs(X) sqrt((X).y*(X).y+(X).x*(X).x)
#endif
#define CMPLX(X,real,imag) {(X).x=(real);(X).y=(imag);}
#define CASSN(to,from) {to.x=from->x;to.y=from->y;}
#define CTREAL(z,X,real) {(z).x=(X).x*(real);(z).y=(X).y*(real);}
#define CSET(to,from) {to->x=(from).x;to->y=(from).y;}