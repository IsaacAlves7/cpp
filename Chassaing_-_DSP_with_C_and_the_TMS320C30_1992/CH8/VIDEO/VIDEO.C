/*VIDEO.C-PROGRAM FOR IMAGE PROCESSING (VIDEO LINE RATE) ANALYSIS*/
#include <c:\c30hll\math.h>
#include <c:\c30hll\stdlib.h>
#include "c30_1.h"
#include "projcmd.h"
#include "serial1.c"
#define buffer 512
#define STAGES 3         /* number of 2nd-order stages in filter */
int buffer1 = buffer/2;
int dma_int2 = 0x00040000;
volatile int *bus = (volatile int *) 0x808060;
volatile int *dma = (volatile int *) 0x808000;
volatile int *host = (volatile int *) 0x804000;
int index = 0, loop, ramp = 0, arrow = 0;
int *input, *intermediate, *average;

float A0[STAGES][3] = { /* numerator coefficients */
                           {  2.186845E-02, 2.186845E-02, 0.000000E+00 },
                           {  2.186845E-02, 4.373690E-02, 2.186845E-02 },
                           {  2.186845E-02, 4.373690E-02, 2.186845E-02 }};
float B0[STAGES][2] = { /* denominator coefficients */
                           { -8.047146E-01, 0.000000E+00 },
                           { -1.665480E+00, 7.049447E-01 },
                           { -1.832568E+00, 8.759922E-01 }};

float A1[STAGES][3] = { /* numerator coefficients */
                           {  6.308639E-02, 6.308639E-02, 0.000000E+00 },
                           {  6.308639E-02, 1.261727E-01, 6.308639E-02 },
                           {  6.308639E-02, 1.261727E-01, 6.308639E-02 }};
float B1[STAGES][2] = { /* denominator coefficients */
                           { -6.401314E-01, 0.000000E+00 },
                           { -1.356731E+00, 4.939717E-01 },
                           { -1.608209E+00, 7.708878E-01 }};

float A2[STAGES][3] = { /* numerator coefficients */
                           {  2.989367E-01, 2.989367E-01, 0.000000E+00 },
                           {  2.989367E-01, 5.978735E-01, 2.989367E-01 },
                           {  2.989367E-01, 5.978735E-01, 2.989367E-01 }};
float B2[STAGES][2] = { /* denominator coefficients */
                           { -1.256041E-01, 0.000000E+00 },
                           { -2.772673E-01, 1.211474E-01 },
                           { -3.806422E-01, 5.391505E-01 }};
float A[STAGES][3];
float B[STAGES][2];
float DLY[STAGES][2] = {0};  /* delay storage elements */
/* 6th order IIR filter */
void IIR(int *IO_in, int *IO_out, int n, int len)
{
  int i, loop = 0, loop1, int_val;
  float dly, yn, input;
  while (loop < len)
  {
    ++loop;
    input = IO_in[loop] & 0xFF;
    for (i = 0; i < n; i++)
    {
      dly = input - B[i][0] * DLY[i][0] - B[i][1] * DLY[i][1];
      yn = A[i][2] * DLY[i][1] + A[i][1] * DLY[i][0] + A[i][0] * dly;
      DLY[i][1] = DLY[i][0];
      DLY[i][0] = dly;
      input = yn;
    }
    if (yn >= 0) int_val = (int)yn;
    else int_val = 0;                       /* clip any negative values */                                                               
    IO_out[loop] = int_val;
  }
  DLY[0][0] = 0;
  DLY[0][1] = 0;
  DLY[1][0] = 0;
  DLY[1][1] = 0;
  DLY[2][0] = 0;
  DLY[2][1] = 0;
}
/* function to perform averaging */
void AVG(int *IO_in, int *IO_out, int loop, int rate)
{
  int index, temp;
  if (loop == 0)
    for (index = 0; index < buffer; index++)
      IO_out[index] = IO_in[index];
  else 
    for (index = 0; index < buffer; index++)
      IO_out[index] += IO_in[index];
  if((rate-loop) == 1)
    for (index = 0; index < buffer; index++) 
      IO_out[index] = IO_out[index] / rate;
} 

main()
{
  int i;
  int *temp;
  int avg = 0;
  int lpa = 0;
  int lpb = 0;
  int lpc = 0;
  int edge1 = 0;
  int edge2 = 0;
  int RATE = 10;
  int line = 128;
  int sel_rate = 1;
  int status;
  int tmp;
  input = (int *) calloc(buffer, sizeof(int));
  average = (int *) calloc(buffer, sizeof(int));
  intermediate = (int *) calloc(buffer1, sizeof(int));
  init_evm();
  asm("  LDI  400h,IE");
  asm("  OR  2000h,ST");
  SERIAL1SET();
  PBASE[0x58] = line;                  /* select line number to be sampled */
  for (;;)
  { 
    for(loop = 0; loop < RATE; loop++) /* input data from SP1 */
    {
      for (index = 0; index < buffer1; index++)
      {
        RWAIT();
        input[index] = PBASE[0x5C];
      }
      for (index = buffer1-1; index >= 0; index--)  /* amplify and shift data*/
      {
        tmp = index*2;
        input[tmp+1] = ((input[index] & 0xFF) - 64) * 1.25;
        input[tmp] = (((input[index] >> 8) & 0xFF) - 64) * 1.25;
      }    
      if (avg)                        /* call averaging function */
        AVG((int *)input, (int *)average, loop, RATE);
      else if (lpa)                    /* select 500 KHz low pass filter */
      {
        for (index = 0; index < 3; index++)
        {
          A[index][0] = A0[index][0];
          B[index][0] = B0[index][0];
        }
        for (index = 0; index < 3; index++)
        {
          A[index][1] = A0[index][1];
          B[index][1] = B0[index][1];
        }
        for (index = 0; index < 3; index++)
        {
          A[index][2] = A0[index][2];
        }
        IIR((int *)input, (int *)average, STAGES, buffer);
      }
      else if (lpb)                    /* select 1 MHz low pass filter */
      {
        for (index = 0; index < 3; index++)
        {
          A[index][0] = A1[index][0];
          B[index][0] = B1[index][0];
        }
        for (index = 0; index < 3; index++)
        {
          A[index][1] = A1[index][1];
          B[index][1] = B1[index][1];
        }
        for (index = 0; index < 3; index++)
        {
          A[index][2] = A1[index][2];
        }
        IIR((int *)input, (int *)average, STAGES, buffer);
      }
      else if (lpc)                    /* select 3 MHz low pass filter */
      {
        for (index = 0; index < 3; index++)
        {
          A[index][0] = A2[index][0];
          B[index][0] = B2[index][0];
        }
        for (index = 0; index < 3; index++)
        {
          A[index][1] = A2[index][1];
          B[index][1] = B2[index][1];
        }
        for (index = 0; index < 3; index++)
        {
          A[index][2] = A2[index][2];
        }
        IIR((int *)input, (int *)average, STAGES, buffer);
      }
      else if (edge1)                  /* first edge enhancement algorithm */
        for (index = 1; index < buffer-1; index++)
        {
          average[index] = 3*input[index]-input[index-1]-input[index+1];  
          if (average[index] > 255) average[index] = 255;
          if (average[index] < 0) average[index] = 0;
        } 
      else if (edge2)                  /* second edge enhancement algorithm */
        for (index = 1; index < buffer-1; index++)
        {
          average[index] = (4*input[index]-input[index-1]-input[index+1])/2;  
          if (average[index] > 255) average[index] = 255;
          if (average[index] < 0) average[index] = 0;
        }
    }
    /* clip data to 8 bits */
    for (index = 0; index < buffer; index++)
    {
      if (input[index] < 0) input[index] = 0;
      if (input[index] > 255) input[index] = 255;
      if (average[index] < 0) average[index] = 0;
      if (average[index] > 255) average[index] = 255;
    }
    /* prepare to send data to PC HOST with 16 bit transfer */
    for (index = 0; index < buffer1; index++) 
    {
      tmp = index*2;
      input[index] = input[tmp] << 8;
      average[index] = average[tmp] << 8;
      input[index] |= input[tmp+1] & 0xFF;
      average[index] |= average[tmp+1] & 0xFF;
    }    
    PBASE[0x58] = line;     /* select next line number to be sampled */
    status = 0;             /* status of EVM to be sent to PC HOST*/  
    if (lpa) status += 1;   /* 500 KHz low pass filter status */
    if (lpb) status += 2;   /* 1 MHz low pass filter status */
    if (lpc) status += 4;   /* 3 MHz low pass filter status */
    if (avg) status += 8;   /* averaging status */
    if (edge1) status += 16;  /* first edge enhancement status */
    if (edge2) status += 32;  /* second edge ehnancement status */
    input[buffer1-2] = status << 8;
    input[buffer1-2] |= line & 0xFF;  /* number of the line sampled */
    input[buffer1-1] = RATE << 8;     /* current rate of the display */
    input[buffer1-1] |= sel_rate & 0xFF;  /* # lines to skip for line select*/
    average[buffer1-2] = input[buffer1-2]; 
    average[buffer1-1] = input[buffer1-1];
    while(dma[TRANSFER]); 
    if (!avg && !lpa && !lpb && !lpc && !edge1 && !edge2)
      temp = input;
    else temp = average;
    intermediate = temp;
    i = *host & 0xFF;
    if (i && i < SEND)
    {
      *host = i;
      switch(i)
      {
        case RESET: {avg = 0; lpa = 0; RATE = 10; arrow = 0; lpb = 0;
                     lpc = 0; sel_rate = 1; edge1 = 0; edge2 = 0; break;}
        case INCREASE: if (RATE < 120) {RATE = RATE << 1; break;}
                       else break; 
        case DECREASE: if (RATE > 1) {RATE = RATE >> 1; break;}
                       else break;
        case LPA: if (lpa) {lpa=0;avg=0;lpb=0;lpc=0;edge1=0;edge2=0;break;}
                  else {lpa=1;avg=0;lpb=0;lpc=0;edge1=0;edge2=0;break;} 
        case LPB: if (lpb) {lpb=0;avg=0;lpa=0;lpc=0;edge1=0;edge2=0;break;}
                  else {lpb=1;avg=0;lpa=0;lpc=0;edge1=0;edge2=0;break;}
        case LPC: if (lpc) {lpc=0;avg=0;lpa=0;lpb=0;edge1=0;edge2=0;break;}
                  else {lpc=1;avg=0;lpa=0;lpb=0;edge1=0;edge2=0;break;}
        case AVER: if (avg) {avg=0;lpa=0;lpb=0;lpc=0;edge1=0;edge2=0;break;} 
                   else {avg=1;lpa=0;lpb=0;lpc=0;edge1=0;edge2=0;break;}
        case EDGE1: if (edge1) {lpc=0;avg=0;lpa=0;lpb=0;edge1=0;break;}
                    else {edge2=0;edge1=1;lpc=0;avg=0;lpa=0;lpb=0;break;}
        case EDGE2: if (edge2) {lpc=0;avg=0;lpa=0;lpb=0;edge2=0;break;}
                    else {edge1=0;edge2=1;lpc=0;avg=0;lpa=0;lpb=0;break;}
        case L_MINUS: if (line > (0 + sel_rate)) {line -= sel_rate; break;}
                      else break;
        case L_PLUS: if (line < (255 - sel_rate)) {line += sel_rate; break;}
                     else break;
        case SEL_RATE: if (sel_rate < 16) {sel_rate = sel_rate << 1; break;}
                       else {sel_rate = 1; break;}
      }
      if (RATE == 4) RATE = 5;
      if (RATE == 8) RATE = 10;
      if (RATE == 15) RATE = 20; 
      if (RATE == 40) RATE = 30;
      while (*host & 0xFF);
      *host = NONE;
    } 
    else if (i == SEND)
    {
      for (loop = 0;loop < buffer;loop++) /*discard field for synchronization*/  
      {                                        
        RWAIT();
        temp[0] = PBASE[0x5C];
      }
      recording((int *)intermediate, SEND); /* transfer results to PC HOST */
    }
  }
}  
    
void init_evm()
{
  bus[EXPANSION] = 0x0;
  bus[PRIMARY] = 0x1000;
  *host = NONE;
  asm("  OR      800h,ST");
}

void c_int11()
{
  while(*host & 0xFF);
  *host = NONE;
  dma[GLOBAL] = 0;
  asm("  AND       0FFFFh,IE");
}

void recording(int *source, int rec_cmd)
{
  asm("  AND     0FFF8h,IF");
  asm("  OR      @_dma_int2,IE");
  dma[SOURCE] = (int) source;
  dma[DEST] = (int) host;
  dma[TRANSFER] = buffer1;
  dma[GLOBAL] = 0x0E13;
  *host = rec_cmd;
}

void c_int05()
{
  for (;;);
}

void c_int99()
{
  for (;;);
}  




