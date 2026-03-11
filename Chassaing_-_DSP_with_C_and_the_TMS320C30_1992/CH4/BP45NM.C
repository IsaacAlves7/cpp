/*BP45NM.C-FIR BANDPASS FILTER IN C WITHOUT USING MODULO*/
#include "bp45coef.h"   /*include coefficient file      */
float DLY[N];           /*delay samples                 */
void filt(float *, float *, int *, int *, int);

main ()
{                    
  int i;
  volatile int *IO_INPUT = (volatile int *) 0x804000;
  volatile int *IO_OUTPUT = (volatile int *) 0x804001;
  for (i = 0; i<N; i++) DLY[i] = 0.0;
  filt((float *)H, (float *)DLY, (int *)IO_INPUT, (int *)IO_OUTPUT, N);
}

void filt(float *h, float *dly, int *IO_input, int *IO_output, int n)
{
  int i, j, m, N1, N1_m, n_m;
  float acc = 0;
  N1 = n-1;
  dly[0] = *IO_input;
  for (m = 0; m < n; m++)
  {        
    N1_m = N1-m;
    n_m = n-m;
    for (i = 0; i < n_m; i++)
      acc += h[i] * dly[N1_m-i];
    for (j = m; j > 0; j--)
      acc += h[n-j] * dly[N1_m+j];
    *IO_output = acc;
    acc = 0.0;
    dly[N1_m] = *IO_input;
  }
}    


