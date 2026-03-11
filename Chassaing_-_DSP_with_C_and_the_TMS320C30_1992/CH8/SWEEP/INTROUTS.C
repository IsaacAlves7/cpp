/* introuts.c */
/* This file contains the interrupt routines for */
/* interrupts 5, 11, and 99 (spurious interrupt) */
                          
#define NONE 0x0
#define GLOBAL 0x0

extern volatile int *host;
extern volatile int *dma;

void c_int05()
{
  for(;;);
}

void c_int11()
{
  while(*host & 0xFF);
  *host = NONE;
  dma[GLOBAL] = 0;

  asm("  AND  0FFFFh,IE");
}

void c_int99()
{
  for (;;);
}

