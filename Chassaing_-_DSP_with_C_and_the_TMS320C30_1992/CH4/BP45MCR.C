/*BP45MCR.C-REAL-TIME FIR BANDPASS FILTER.CALLS ASSEMBLY FUNCTION  */
#include "aiccom.c"                   /*include AIC com routines   */
#include "bp45coef.h"                 /*include coefficients file  */
#define VEC_ADDR (volatile int *) 0x00
float DLY[2*N];                       /*delay samples              */
int AICSEC[4] =  {0x1428,0x1,0x4A96,0x67}; /*AIC config data       */
int data_in, data_out;
extern void filt(float *, float *, int *, int *, int);

void c_int05()
{
  PBASE[0x48] = data_out << 2;
  data_in = PBASE[0x4C] << 16 >> 18;
}

main ()
{
  volatile int *INTVEC = VEC_ADDR;
  int *IO_INPUT, *IO_OUTPUT;
  IO_INPUT = &data_in;
  IO_OUTPUT = &data_out;
  INTVEC[5] = (volatile int) c_int05;
  AICSET_I();
  for (;;) 
    filt((float *)H, (float *)DLY, (int *)IO_INPUT, (int *)IO_OUTPUT, N);
}


