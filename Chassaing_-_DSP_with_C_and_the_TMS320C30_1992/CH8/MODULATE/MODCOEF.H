/*MODCOEF.H COEFFICIENTS FOR 1,2,3,4 kHz BANDPASS FILTERS*/
#define NUM_STAGES 3         /* number of 2nd-order stages in filter */
double A1k[NUM_STAGES][3] = { /* numerator coefficients */
       {  1.06900169486E-0001, -1.80934638738E-0001,  1.06900169486E-0001 },
       {  1.06900169486E-0001, -1.63523188763E-0001,  1.06900169486E-0001 },
       {  1.06900169486E-0001,  0.00000000000E+0000, -1.06900169486E-0001 }
                          };
double B1k[NUM_STAGES][2] = { /* denominator coefficients */
       { -1.64606071835E+0000,  9.94221858963E-0001 },
       { -1.60715196747E+0000,  9.85568839917E-0001 },
       { -1.57928976080E+0000,  9.93754126054E-0001 }
                          };
double DLY1k[NUM_STAGES][2];  /* delay storage elements (z-shifts) */
double A2k[NUM_STAGES][3] = { /* numerator coefficients */
       {  1.05031366006E-0001, -7.87438109303E-0002,  1.05031366006E-0001 },
       {  1.05031366006E-0001, -5.04973319662E-0002,  1.05031366006E-0001 },
       {  1.05031366006E-0001,  0.00000000000E+0000, -1.05031366006E-0001 }
                          };
double B2k[NUM_STAGES][2] = { /* denominator coefficients */
       { -6.69667282668E-0001,  9.94009116616E-0001 },
       { -6.13895264178E-0001,  9.85626535409E-0001 },
       { -5.62287536701E-0001,  9.93898189038E-0001 }
                          };
double DLY2k[NUM_STAGES][2];  /* delay storage elements (z-shifts) */
double A3k[NUM_STAGES][3] = { /* numerator coefficients */
       {  1.05031366006E-0001,  7.87438109303E-0002,  1.05031366006E-0001 },
       {  1.05031366006E-0001,  5.04973319662E-0002,  1.05031366006E-0001 },
       {  1.05031366006E-0001,  0.00000000000E+0000, -1.05031366006E-0001 }
                          };
double B3k[NUM_STAGES][2] = { /* denominator coefficients */
       {  6.69667282668E-0001,  9.94009116616E-0001 },
       {  6.13895264178E-0001,  9.85626535409E-0001 },
       {  5.62287536701E-0001,  9.93898189038E-0001 }
                          };
double DLY3k[NUM_STAGES][2];  /* delay storage elements (z-shifts) */
double A4k[NUM_STAGES][3] = { /* numerator coefficients */
       {  1.06900169486E-0001,  1.80934638738E-0001,  1.06900169486E-0001 },
       {  1.06900169486E-0001,  1.63523188763E-0001,  1.06900169486E-0001 },
       {  1.06900169486E-0001,  0.00000000000E+0000, -1.06900169486E-0001 }
                          };
double B4k[NUM_STAGES][2] = { /* denominator coefficients */
       {  1.64606071835E+0000,  9.94221858963E-0001 },
       {  1.60715196747E+0000,  9.85568839917E-0001 },
       {  1.57928976080E+0000,  9.93754126054E-0001 }
                          };
double DLY4k[NUM_STAGES][2];  /* delay storage elements (z-shifts) */





