/*BP45ERIC.C-FIR BANDPASS WITH 45 COEFFICIENTS. SAMPLES SHIFTED*/
#include "bp45coef.h"           /*include coefficient file     */
#define N 45                    /*length of impulse response   */
float DLY[N*2-1];               /*init for 2*N-1 samples       */
void filt(float *,float *,int *,int *,int); /*filter routine   */ 

main()
{
  volatile int *in_addr =(volatile int *)0x804000; /*input addr */
  volatile int *out_addr=(volatile int *)0x804001; /*output addr*/       
  filt((float *)H,(float *)DLY,(int *)in_addr,(int *)out_addr,N);
}                                                               

void filt(float *h,float *dly,int *in_port,int *out_port,int NS)       
{
  float acc=0.0;             /*init accumulator         */
  int   i,j,k=NS-1;          /*index variables          */ 
  for (i = 0; i<N*2-1; i++) DLY[i] = 0.0; /*init samples*/
  for (i=0;i<NS;i++)
  {
     dly[i+k] =*in_port;     /*get new sample           */
     for (j=0;j<NS;j++)                        
       acc += h[j]*dly[i+j]; /*perform convolution      */       
     *out_port=acc;          /*output new value         */
     acc=0.0;                /*reset accumulator        */
  }                      
  for (i=0;i<k;i++)          /*shift values from        */
    dly[i]=dly[i+NS];        /*lower half to upper half */
}






























