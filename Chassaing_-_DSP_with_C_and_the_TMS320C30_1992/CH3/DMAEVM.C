/*DMAEVM.C-HOST/EVM COMMUNICATION WITH DMA*/
#include "c30_1.h"
#define buffer 10
#define SEND 128
int dma_int1 = 0x00020000;
int dma_int2 = 0x00040000;
volatile int *bus = (volatile int *) 0x808060;
volatile int *dma = (volatile int *) 0x808000;
volatile int *host = (volatile int *) 0x804000;

void main(void)
{
  int i;
  int data[10];
  int temp[10];
  asm("  LDI  400h,IE");
  asm("  LDI  2000h,ST");
  *host = 0x00;
  do                        
  {
    i = *host & 0xFF;            
    if (i == 64)
    {
      playing (data, 64);    /*transfer 10 values from host*/
      while(dma[TRANSFER]);
    }
  }
  while (i != 64);
  for (i = 0; i < 10; i++)
    temp[i] = data[9-i];     /*reverse order of values*/
  do
  {
    i = *host & 0xFF;
    if (i == SEND)
      recording(temp, SEND); /*transmit values to host*/ 
  }                          /*in reversed order*/
  while (i != 128); 
}  

void recording(int *source, int rec_cmd)
{
  asm("  AND    0FFF8h,IF");
  asm("  OR     @_dma_int2,IE");
  dma[SOURCE] = (int)source;
  dma[DEST] = (int)host;
  dma[TRANSFER] = buffer;
  dma[GLOBAL] = 0x0E13;
  *host = rec_cmd;
}

void playing(int *dest, int play_cmd)
{
  asm("  AND    0FFF8h,IF");
  asm("  OR     @_dma_int1,IE");
  dma[SOURCE] = (int)host;
  dma[DEST] = (int)dest;
  dma[TRANSFER] = buffer;
  dma[GLOBAL] = 0x0D43;
  *host = play_cmd;
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






