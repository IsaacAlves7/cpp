/*BP45NMR.C-REAL-TIME FIR BANDPASS FILTER IN C.CALLS AICCOM */
#include "aiccom.c"              /*include AIC com routines */
#include "bp45coef.h"            /*include coefficients file*/
#define VEC_ADDR (volatile int *) 0x00
float DLY[N];                    /*delay samples            */
int data_in, data_out;
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67};  /*AIC config data*/

void filt(float *h, float *dly, int *IO_input, int *IO_output, int n)
{
  int i, j, m, N1, N1_m, n_m, index = 0;
  float acc = 0;
  N1 = n-1;
  dly[0] = *IO_input;
  for (m = 0; m < n; m++)
  {
    asm("     idle");             /*wait for interrupt        */
    N1_m = N1-m;
    n_m = n-m;
    for (i = 0; i < n_m; i++)     /*addr below new sample to 0*/
      acc += h[i] * dly[N1_m-i];
    for (j = m; j > 0; j--)       /*from n to latest sample   */
      acc += h[n-j] * dly[N1_m+j]; /*latest sample last       */ 
    *IO_output = acc;             /*output result             */
    acc = 0.0;                    /*clear accumulator         */
    dly[N1_m] = *IO_input;        /*get new sample            */
  }
}    

void c_int05()
{
  data_in = UPDATE_SAMPLE(data_out);
}

main ()
{
  volatile int *INTVEC = VEC_ADDR;
  int *IO_INPUT, *IO_OUTPUT;
  IO_INPUT = &data_in;
  IO_OUTPUT = &data_out;
  INTVEC[5] = (volatile int) c_int05;
  AICSET_I();
  for(;;)
   filt((float *)H, (float *)DLY, (int *)IO_INPUT, (int *)IO_OUTPUT, N);
}




