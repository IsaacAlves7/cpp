/*BP45H.C-FIR BANDPASS FILTER GENERATED USING HYPERSIGNAL MODIFIED FOR I/O
**************************************************************************
* Texas Inst. TMS320C30 C Lang. Filter Realization Copyright (C) by HYPERCEPTION, INC.
**************************************************************************
*
*       Filter Generated from File: bp45.FIR
*
*       Filter Order=45
*
*       Direct Form Realization of the following convolution sum:
*
*                               N-1
*                               ____
*                               \
*                      y(n) =   /___  h(k)x(n-k)
*                              k=0
*
*
*               -1     -1                 -1
*              z      z                  z
*  o--->---o--->---o--->---o---> - - -o--->----o
*   x(n)   |       |       |          |        |
*          |       |       |          |        |
*          v h(0)  v h(1)  v h(2)     v h(N-2) v h(N-1)
*          |       |       |          |        |
*          |       |       |          |        |
*          o--->-- + -->-- + --> - -  + -->--- + -->---o
*                                                  y(n)
*
*
*/
#define N 45         /* length of FIR filter impulse response */
int start_index = { 0 }; /* circular buffer starting position */
double H[N] = { /* filter coefficients */
 -1.83989092791E-0003, -2.65748927587E-0003, -4.31240164415E-0010, 
  3.15489863431E-0003,  2.59504949565E-0003, -4.15992854430E-0003, 
 -1.54063990002E-0002, -2.50783676799E-0002, -2.54745141379E-0002, 
 -1.17986263157E-0002,  1.39224582794E-0002,  4.20645113719E-0002, 
  5.88840954593E-0002,  5.30715079282E-0002,  2.22550923780E-0002, 
 -2.41095421284E-0002, -6.75430291924E-0002, -8.83147433676E-0002, 
 -7.47571889853E-0002, -2.95655400338E-0002,  3.03027241599E-0002, 
  8.05076375048E-0002,  1.00000000600E-0001,  8.05076375048E-0002, 
  3.03027241599E-0002, -2.95655400338E-0002, -7.47571889853E-0002, 
 -8.83147433676E-0002, -6.75430291924E-0002, -2.41095421284E-0002, 
  2.22550923780E-0002,  5.30715079282E-0002,  5.88840954593E-0002, 
  4.20645113719E-0002,  1.39224582794E-0002, -1.17986263157E-0002, 
 -2.54745141379E-0002, -2.50783676799E-0002, -1.54063990002E-0002, 
 -4.15992854430E-0003,  2.59504949565E-0003,  3.15489863431E-0003, 
 -4.31240164415E-0010, -2.65748927587E-0003, -1.83989092791E-0003};
double DLY[N];  /* delay taps (z-shifts) */
double filt(stage_input)  /* actual filter routine */
double stage_input;
{
    double acc;
    int   i, j;
    DLY[start_index] = stage_input;
    j = --start_index;
    acc = 0.0;
    for (i=0; i<N; i++)
    {
        j = ++j % N;     /* circular buf requires modulo inc */
        acc += H[i] * DLY[j];
    }
    start_index = j;
    return acc;
}

main ()
{
#define IMPULSE_LENGTH 45   /* length of FIR filter impulse response */
volatile int *IO_OUT = (volatile int *) 0x804001; /*added for i/o    */
/* double out_val[IMPULSE_LENGTH];     can be deleted                */
int   n;
/*     out_val[0]=filt(10000.0);  /* the "impulse" input --> deleted */
    *IO_OUT = filt(10000.0);      /* the "impulse"                   */
    for (n=1; n<IMPULSE_LENGTH; n++)
    {
     *IO_OUT = filt(0.0);         /* other values are zero           */
    }
}


