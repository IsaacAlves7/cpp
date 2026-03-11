/*FFT8MC.C-8-POINT FFT.CALLS FUNCTION FFT_RL.ASM IN TMS320C30 CODE*/
#include <math.h>
#define N 8        /*FFT length  */
#define M 3        /*# of stages */
float data[N] = {1,1,1,1,0,0,0,0};   /*real-valued input samples*/   
float real1, img1;
extern void fft_rl(int, int, float *);  /*generic FFT function  */
volatile int *IO_OUT = (volatile *) 0x804001; /*output port addr*/  

main()
{
  int loop;
  fft_rl(N, M, (float *)data); 
  *IO_OUT = (int)(data[0]*1000);
  for (loop = 1; loop < N/2; loop++)
  {
    real1 = data[loop];
    img1 = data[N-loop];
    *IO_OUT = (int)(real1*1000);
    *IO_OUT = (int)(img1*1000);
  }
  *IO_OUT = (int)(data[N/2]*1000);
  for (loop = N/2+1; loop < N; loop++)
  {
    real1 = data[N-loop];
    img1 = data[loop];
    *IO_OUT = (int)(real1*1000);
    *IO_OUT = (int)(img1*(-1000));
  }
}



