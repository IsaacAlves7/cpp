/*FFT512R.C-REAL-TIME FFT WITH 512 POINTS. CALLS FFT_RL*/
#include <math.h>                  /*std library func  */
#include "aiccom.c"                /*AIC com routines  */
#define N 512                      /*size of FFT       */
#define M 9                        /*number of stages  */
volatile int index = 0;            /*input_output index*/
float *IO_buffer, *data, *temp;          /*pointers to array buffers*/
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67};     /*AIC config data     */
extern void fft_rl(int, int, float *);        /*fft function protype*/
volatile int *INTVEC = (volatile int *) 0x00; /*addr of inter vecs  */
void c_int05()                         /*interrupt handler function */
{
 PBASE[0x48] = ((int)(IO_buffer[index])) << 2;        /*output data */
 IO_buffer[index] = (float)(PBASE[0x4C] << 16 >> 18); /*input data  */
  if (++index >= N) index = 0;         /*increment index, reset = N */
}
main()
{
 int loop;                                       /* declare variable       */
 float real, img;                                /* declare variables      */
 INTVEC[5] = (volatile int) c_int05;             /*install inter function  */
 AICSET_I();                                     /*config AIC for interrupt*/
 IO_buffer = (float *) calloc(N, sizeof(float)); /*input_output buffer     */
 data = (float *) calloc(N, sizeof(float));      /* fft data buffer        */
 while (1)                                       /* create endless loop    */
 {
  fft_rl(N, M, (float *)data);             /*call FFT function       */
  data[0] = sqrt(data[0]*data[0])/N;       /*magnitude of X(0)       */
  for (loop = 1; loop < N/2; loop++)       /*calculate X(1)..X(N/2-1)*/
  {
   real = data[loop];                      /*real part               */
   img = data[N-loop];                     /*imaginary part          */
   data[loop] = sqrt(real*real+img*img)/N; /*find magnitude          */
  }
  data[N/2] = sqrt(data[N/2]*data[N/2])/N; /*magnitude of X(N/2)     */
  for (loop = N/2+1; loop < N; loop++)     /*X(N/2+1).. X(N-1)       */
   data[loop] = data[N-loop];              /*use symetry             */ 
  while (index);                           /*wait till IO_buffer is empty*/
  temp = data;                             /*temp pointer => data buffer */
  data = IO_buffer;                  /*IO_buffer becomes new data buffer */
  IO_buffer = temp;                  /*data buffer becomes new IO_buffer */
  IO_buffer[0] = -2048;              /*frame sync pulse (negative spike) */
 }
}



