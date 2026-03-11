/*BP45NMD.C-FIR BANDPASS FILTER IN C WITHOUT USING MODULO*/
#include "bp45coef.h"    /*include coefficient file      */
float DLY[N];            /*delay samples                 */
void filt(float *, float *, int *, int *, int);

main ()
{
  int i;
  volatile int *IO_INPUT = (volatile int *) 0x804000;
  volatile int *IO_OUTPUT = (volatile int *) 0x804001;
  for (i = 0; i < N; i++) DLY[i] = 0.0;
  filt((float *)H, (float *)DLY, (int *)IO_INPUT, (int *)IO_OUTPUT, N);
}

void filt(float *h, float *dly, int *IO_input, int *IO_output, int n)
{
  int i, t;
  float acc;
  for (t = 0; t < n; t++)
  { 
    acc = 0.0;
    dly[0] = *IO_input;
    for (i = 0; i < n; i++)
      acc += h[i] * dly[i];
    for (i = n-1; i > 0; i--)
      dly[i] = dly[i-1];   
    *IO_output = acc;
  }
}    


