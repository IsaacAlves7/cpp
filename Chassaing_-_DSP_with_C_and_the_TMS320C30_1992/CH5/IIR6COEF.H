/*IIR6COEF.H-COEFF FOR SIXTH-ORDER IIR BANDPASS FILTER  */
#define STAGES 3            /*number of 2nd-order stages*/
float A[STAGES][3]=        {/*numerator coefficients    */
{5.3324E-02, 0.0000E+00, -5.3324E-02},
{5.3324E-02, 0.0000E+00, -5.3324E-02},
{5.3324E-02, 0.0000E+00, -5.3324E-02} };
float B[STAGES][2]=        {/*denominator coefficients  */
{-1.4435E+00, 9.4879E-01},
{-1.3427E+00, 8.9514E-01},
{-1.3082E+00, 9.4377E-01} };
   


