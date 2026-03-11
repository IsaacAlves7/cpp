/*FFT8.C-MAIN FFT.RECTANGULAR INPUT DATA.CALLS FFT.C*/
#include <complex.h>  /*complex structure definition*/
extern void FFT();    /*FFT function                */
volatile int *out_addr=(volatile int *)0x804001;
       
main()
  {             /*complex input samples */
  COMPLEX y[8]={10.0,0.0,10.0,0.0,10.0,0.0,10.0,0.0,
                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  int i, n=8;       

  FFT(y,n);          /*calls generic FFT function*/
  for (i=0;i<n;i++)       
    {
    *out_addr=(y[i]).real;
    *out_addr=(y[i]).imag;
    }    
  }



















