/*IIRA.C - IIR FILTER*/
#define STAGES 3          /* number of 2nd-order stages in filter */
#include <stdio.h>
double A[STAGES][3] = { /* numerator coefficients */
                      {   5.3324E-02,  0.0000E+00, -5.3324E-02},
                      {   5.3324E-02,  0.0000E+00, -5.3324E-02},
                      {   5.3324E-02,  0.0000E+00, -5.3324E-02}};
double B[STAGES][2] = { /* denominator coefficients */
		      {  -1.4435E+00,  9.4879E-01},
                      {  -1.3427E+00,  8.9514E-01},
                      {  -1.3082E+00,  9.4377E-01}};
double DLY[STAGES][2] = {0};        /*delay storage elements*/
FILE *IO_IN, *IO_OUT;

void IIR(int n, int len)            /*IIR function*/
{
  int i, loop = 0;
  double dly, yn, input;
  long temp;
  while (loop < len)                /*while loop*/
  {
    ++loop;
    fscanf(IO_IN, "%D\n", &temp);
    input = temp;
    for (i = 0; i < n; i++)         /*for loop*/
    {
      dly = input - B[i][0]*DLY[i][0] - B[i][1]*DLY[i][1];
      yn = A[i][2]*DLY[i][1] + A[i][1]*DLY[i][0] + A[i][0]*dly;
      DLY[i][1] = DLY[i][0];
      DLY[i][0] = dly;
      input = yn;
    }
    fprintf (IO_OUT, "%d\n", (int)yn);   /*output data to file*/
  }
}

main()                                   /*main function*/
{
  #define length 345
  IO_IN = fopen ("IMPULSE.DAT", "rt");   /*open IMPULSE data file*/
  IO_OUT = fopen ("IIRA.DAT", "wt");     /*open OUTPUT data file */
  IIR(STAGES, length);                   /*execute IIR function  */
  fclose (IO_IN);                        /*close IMPULSE.DAT file*/
  fclose (IO_OUT);                       /*close IIROUT.DAT file */
}

