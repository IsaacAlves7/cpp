/*MOD.C Modulator.Requires 1 kHz carrier input.interactive with MODHOST.C*/
#include "math.h"
#include "aiccom.c"              /*include AIC comm routines*/
#include "modcoef.h"             /*coefficients file        */
#include "c30_1.h"               /*DMA/PC Comm routines    */
#define VEC_ADDR (volatile int *) 0x00
#define PI 3.141592654
#define sampling_freq 10000.0
#define buffer 1                           /*1 word DMA transfer*/
int dma_int1 = 0x00020000;
int dma_int2 = 0x00040000;
volatile int *dma = (volatile int *) 0x808000;
volatile int *host = (volatile int *) 0x804000;
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67};          /*AIC config data*/
int data_in, data_out;
int data[1] = {1000};
int k;
float freq = 1000.0;               /*initialize variable */

void GET_FREQ()
{
  k = *host & 0xFF;            
  if (k == 64)                     /*any read DMA data???  */
    {
      playing (data, 64); 
      while(dma[TRANSFER]);        /*wait until buffer = 0 */
      while(k != 64);              /* wait until complete */
      freq = (float)data[0];
    }
}

void playing(int *dest, int play_cmd)   /*Read DMA port data routine*/
{
  asm("  AND    0FFF8h,IF");
  asm("  OR     @_dma_int1,IE");
  dma[SOURCE] = (int)host;
  dma[DEST] = (int)dest;
  dma[TRANSFER] = buffer;
  dma[GLOBAL] = 0x0D43;
  *host = play_cmd;
}
  
void c_int11()              /*DMA interupt handler*/
{
  while(*host & 0xFF);
  *host = NONE;
  dma[GLOBAL] = 0;
  asm("  AND    0FFFFh,IE");
}

void c_int99()
{
  for (;;);
}

float IIR(int *X,int *Y,int n,int l,double A[][3],double B[][2],double DLY[][2])
{
  int i, loop = 0;
  float dly, yn, input;
  double cosine;
  while (loop < l)
  {
    asm("    IDLE   ");         /*wait for update*/
    ++loop;
    input = *X;
    cosine = cos(2.0*PI*loop*(freq/sampling_freq));  /*modulation eq. */
    input = input*cosine; 
    for (i = 0; i < n; i++)
      {
      dly = input - B[i][0] * DLY[i][0] - B[i][1] * DLY[i][1];
      yn = A[i][2] * DLY[i][1] + A[i][1] * DLY[i][0] + A[i][0] * dly;
      DLY[i][1] = DLY[i][0];
      DLY[i][0] = dly;
      input = yn;
      }
    *Y = yn/8.0 ;     /*cancel AIC amplification factor*/
  }
}

void c_int05()  
{
  PBASE[0x48] = data_out << 2;
  data_in = PBASE[0x4C] << 16 >> 18;
}

main()
{
  volatile int *INTVEC = VEC_ADDR;
  int *IO_INPUT, *IO_OUTPUT,length;
  float four_cycles;
  IO_INPUT = &data_in;
  IO_OUTPUT = &data_out;
  INTVEC[5] = (volatile int) c_int05;
  asm("  LDI  400h,IE");
  asm("  LDI  2000h,ST");
  AICSET_I();
  *host = 0x00;
  for (;;) 
  {
   GET_FREQ();      /*check if Host has new shift frequency*/
   length = 1000;
   if(freq==1.0)
   {
    IIR((int *)IO_INPUT,(int *)IO_OUTPUT,NUM_STAGES,length,A1k,B1k,DLY1k);
   }
   if(freq==1000.0)
   {
    IIR((int *)IO_INPUT,(int *)IO_OUTPUT,NUM_STAGES,length,A2k,B2k,DLY2k);
   }
   if(freq==2000.0)
   {
    IIR((int *)IO_INPUT,(int *)IO_OUTPUT,NUM_STAGES,length,A3k,B3k,DLY3k);
   }
   if(freq==3000.0)
   {
    IIR((int *)IO_INPUT,(int *)IO_OUTPUT,NUM_STAGES,length,A4k,B4k,DLY4k);
   }
  }
}        
   


