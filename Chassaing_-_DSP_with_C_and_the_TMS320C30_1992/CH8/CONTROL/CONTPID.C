/*CONTPID.C-PROGRAM FOR PID CONTROLLER*/
#include "aiccom01.c"
#include "c30_1.h"
#define buffer 3
#define SEND 128
int AICSEC[4];
int AICSEC1[4] = {0x2040,0x1,0x62C6,0x63};
int dma_int1 = 0x00020000;
int dma_int2 = 0x00040000;
volatile int *bus = (volatile int *) 0x808060;
volatile int *dma = (volatile int *) 0x808000;
volatile int *host = (volatile int *) 0x804000;

void main(void)
{
  int i;
  int loop = 1;
  int data[3] = {0,0,0};
  float CNTRL;
  float MODEL = 0;
  float OMEGA;
  float PCNTRL = 0.0;
  float ICNTRL = 0.0;
  float DCNTRL = 0.0;
  float GOAL = 1400; 
  float PGAIN = 0.0;
  float IGAIN = 0.0;
  float DGAIN = 0.0;
  float ERROR;
  float ERROR1 = 0.0;
  float OFFSET = 2000;
  init_evm();
  asm("  LDI  400h,IE");
  asm("  OR  2000h,ST");
  AICSET(SP1);
  PBASE[0x38] = 0x00000008; 
  do
  {
    for (i = 0; i < 5; i++)
      OMEGA = UPDATE_SAMPLE(SP1, 0 - OFFSET);
      OMEGA = -1 * OMEGA;
    i = *host & 0xFF;
    if (i == 64)
    {
      playing (data, 64);
      while(dma[TRANSFER]);
    }
  }
  while (i != 64);
  PGAIN = (float)(data[0] & 0xFFFF);
  IGAIN = (float)(data[1] & 0xFFFF)/1000;
  DGAIN = (float)(data[2] & 0xFFFF);
  do
  {
    for (i = 0; i < 20; i++)
    {
      if (MODEL <= GOAL) MODEL += 1;
      ERROR = MODEL - OMEGA;
      DCNTRL = (ERROR - ERROR1) * DGAIN; 
      PCNTRL = ERROR * PGAIN;
      ICNTRL += ERROR * IGAIN;
      CNTRL = PCNTRL + ICNTRL + DCNTRL;
      if (CNTRL > 4000) CNTRL = 4000;
      if (CNTRL < 0) CNTRL = 0;
      ERROR1 = ERROR;
      OMEGA = UPDATE_SAMPLE(SP1, CNTRL - OFFSET) + 60;
      OMEGA = (-1 * OMEGA);
    }
    data[0] = MODEL*(200/GOAL);
    data[1] = OMEGA*(200/GOAL);
    data[2] = CNTRL/64;
    i = *host & 0xFF;
    if (i == SEND)
      recording(data, SEND);   
  }
  while (1);
}  

void c_int11()
{
  while(*host & 0xFF);
  *host = NONE;
  dma[GLOBAL] = 0;
  asm("  AND    0FFFFh,IE");
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

void c_int05()
{
  for (;;);
}

void c_int99()
{
  for (;;);
}

void init_evm()
{
  bus[EXPANSION] = 0x0;
  bus[PRIMARY] = 0x1000;
  *host = NONE;
  asm("  OR   800h,ST");
}






