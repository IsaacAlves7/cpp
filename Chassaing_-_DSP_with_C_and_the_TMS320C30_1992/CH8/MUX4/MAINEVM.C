/*MAINEVM.C-PROGRAM FOR 4-CHANNEL MULTIPLEXER*/
#include "c30_1.h"
#include "math.h"
#define buffer1 1
#define N 128
#define M 7
#define SEND 128
int dma_int1 = 0x00020000;
int dma_int2 = 0x00040000;
int *serglob1     = (int *)0x808050L; /* Serial Port 1 Global Control Register  */
int *serprtx1     = (int *)0x808052L; /* Serial Port 1 FSX/DX/CLKX Port Control */
int *serprtr1     = (int *)0x808053L; /* Serial Port 1 FSR/DR/CLKR Port Control */
int *sertim1      = (int *)0x808054L; /* Serial 1 Timer Control Register */
int *sertimcntr1  = (int *)0x808055L; /* Serial 1 Timer Counter Register */
int *sertimper1   = (int *)0x808056L; /* Serial 1 Timer Period Register */
int *serxmit1     = (int *)0x808058L; /* Serial 1 Data Transmit Register */
int *serrxd1      = (int *)0x80805CL; /* Serial 1 Data Receive Register */
int ad_data[N]; 
int chan[1];
float *IO_buffer,*data;
volatile int *bus = (volatile int *) 0x808060;
volatile int *dma = (volatile int *) 0x808000;
volatile int *host = (volatile int *) 0x804000;
extern void fft_rl(int, int, float *);
void sp_init(int *chan);
void sp_xfer(int *datbuf);
void c_int11();
void c_int05();
void c_int99();

void main(void)
{
  int i,loop;
  float real, img,res;
  IO_buffer = (float *) calloc(N, sizeof(float));
  data = (float *) calloc(N, sizeof(float));
  asm("  LDI  400h,IE");
  asm("  OR   2000h,ST");
 for(;;)
 {
  *host = 0x00;
  do                        
  {
    i = *host & 0xFF;            
    if (i == 64)
    {
      playing (chan, 64);
      while(dma[TRANSFER]);
    }
  }
  while (i != 64);
   sp_init(chan); 
   sp_xfer(ad_data);
   for(i = 0; i < N; i++)
	data[i] = (float)(ad_data[i]);

   fft_rl(N, M, (float *)data);
   
   data[0] = sqrt(data[0]* data[0])/N;
   for (loop = 1; loop < N/2; loop++)      
   {
    real = data[loop];
    img = data[N-loop];
    res = real * real + img * img;
    data[loop] = sqrt(res)/N;
   }
   data[N/2] = sqrt(data[N/2]* data[N/2])/N;
   for (loop = N/2+1; loop < N; loop++)   
	data[loop] = data[N-loop];
   data[0] = 0.0;
  for (i = 0; i < N; i++)
    ad_data[i] = (int)data[i];
  do                            /* send it to the pc */
  {
    i = *host & 0xFF;
    if (i == 128)
      recording((int *)ad_data, 128);
  }                          
  while (i != 128); 
 }
}  

void recording(int *source, int rec_cmd)  /* send data out to pc */
{
  asm("  AND    0FFF8h,IF");
  asm("  OR     @_dma_int2,IE");
  dma[SOURCE] = (int)source;
  dma[DEST] = (int)host;
  dma[TRANSFER] = N ;   /*  N/2 */
  dma[GLOBAL] = 0x0E13;
  *host = rec_cmd;
}

void playing(int *dest, int play_cmd)     /*  recieve data from pc  */
{
  asm("  AND    0FFF8h,IF");
  asm("  OR     @_dma_int1,IE");
  dma[SOURCE] = (int)host;
  dma[DEST] = (int)dest;
  dma[TRANSFER] = buffer1;
  dma[GLOBAL] = 0x0D43;
  *host = play_cmd;
}
void sp_init(int *ch)
{
   int i,n;
   *serprtx1   = 0X00000111;     /* xmit port control register   */
   *serprtr1   = 0X00000111;     /* rec port control register    */
   *sertimper1 = 0x00080008;            /* timer period register */
   *sertim1    = 0x3cf;          /* timer control register was 0x0f;*/
   *serglob1   = 0XC0000C4;      /* was  *serglob1   = 0X4000144;   */
   for (i = 0; i < 10; i++);
   while(!(*serglob1 & 0x02));
   *serxmit1 = 0x00;
   for (i = 0; i < 10; i++);
   while(!(*serglob1 & 0x02));
   *serxmit1 = ch[0];
}
void sp_xfer(int *datbuf)
{
   int i,value = 0;
   int temp_dat[256];
   for (i = 0; i < 256; i++)
   {
	while(!(*serglob1 & 0x01));
	value = (*serrxd1 >> 8) & 255;
	temp_dat[i] = value;
    }
   for (i = 0; i < N; i++)
      datbuf[i] = temp_dat[i+5];
}
  
void c_int11()
{
  while(*host & 0xFF);
  *host = NONE;
  dma[GLOBAL] = 0;
  asm("  AND    0FFFFh,IE");
}

void c_int05()
{
  for (;;);
}

void c_int99()
{
  for (;;);
}







