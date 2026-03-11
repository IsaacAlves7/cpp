#include "Cmatrix.hpp"

#define maxterm 100
#define abstol 1.e-10
#define reltol 1.e-5
#define errorcode -1.e60
#define pi 3.141592653589793238462643383279
#define lim 15.
#define max(a,b) ((a)>(b)?(a):(b))
#define infinite_loop for(;;)

//int incomplete;
complex agm( complex& ,complex& ,complex& );
complex F( complex& ,complex& );
complex K(complex& );
complex E(complex& );
complex Z(complex& ,complex& );
complex E(complex& ,complex& );
extern complex jtheta,amplitude;
void scdn( complex& ,complex& ,complex& ,complex& ,complex& );

complex heuman(complex& phi,complex& m);
void theta(complex& v,complex& m,
complex& t1,complex& t2,complex& t3,complex& t4);

void neville(complex& u,complex& m,
complex& ts,complex& tc,complex& td,complex& tn) ;

complex Weierstrass( complex&, complex&, complex&);
complex InvWeier( complex&, complex&, complex&);
complex ElLog( complex&, complex&, complex&);
complex ElExp( complex&, complex&, complex&);

complex tangent ( complex& x);
complex tanh(complex& x);
complex sinh(complex& x);
complex asin(complex&x);
complex acos(complex&x);
complex arctan(complex&x);
