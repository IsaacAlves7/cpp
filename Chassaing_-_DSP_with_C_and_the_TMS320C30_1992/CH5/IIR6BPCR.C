/*IIR6BPCR.C-REAL-TIME IIR BANDPASS FILTER IN C   */
#include "aiccom.c"    /*include AIC comm routines*/
#include "iir6coef.h"  /*coefficients file        */
#define VEC_ADDR (volatile int *) 0x00

float DLY[STAGES][2] = {0}; /*delay samples       */
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67}; /*AIC config data*/
int data_in, data_out;

float IIR(int *IO_in, int *IO_out, int n, int len)
{
  int i, loop = 0;
  float dly, yn, input;

  while (loop < len)
  {
    asm("    IDLE   ");
    ++loop;
    input = *IO_in;   
    for (i = 0; i < n; i++)
      {
      dly = input - B[i][0] * DLY[i][0] - B[i][1] * DLY[i][1];
      yn = A[i][2] * DLY[i][1] + A[i][1] * DLY[i][0] + A[i][0] * dly;
      DLY[i][1] = DLY[i][0];
      DLY[i][0] = dly;
      input = yn;
      }
    *IO_out = yn;
  }
}

void c_int05()
{
  PBASE[0x48] = data_out << 2;
  data_in = PBASE[0x4C] << 16 >> 18;
}

main()
{
  #define length 345
  volatile int *INTVEC = VEC_ADDR;
  int *IO_OUTPUT, *IO_INPUT;
  IO_INPUT = &data_in;
  IO_OUTPUT = &data_out;
  INTVEC[5] = (volatile int) c_int05;
  AICSET_I();
  for (;;)  
    IIR((int *)IO_INPUT, (int *)IO_OUTPUT, STAGES, length);
}        
   
