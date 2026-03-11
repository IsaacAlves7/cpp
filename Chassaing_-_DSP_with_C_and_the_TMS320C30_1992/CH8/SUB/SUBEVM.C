/*SUBEVM.C TRACKING OF A SUBMERSIBLE.CALLS FFT_RL.ASM*/
#include "math.h"                  /* Standard library functions */
#include "aiccom.c"                /* AIC commmunications routines */
#define N 512                      /* N sets the size of the FFT */
#define M 9                        /* M defines the number of stages */
#include "c30_1.h"
#define buffer 512
#define SEND 128
int dma_int1 = 0x00020000;
int dma_int2 = 0x00040000;
volatile int *dma = (volatile int *) 0x808000;     /* DMA address defined */
volatile int *host = (volatile int *) 0x804000;    /* Host address defined */
volatile int index = 0;                     /* Input_Output index */
float *IO_buffer, *data, *temp;             /* Pointers to array buffers */
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67};   /* AIC configuration data */
int fft_val[N];
extern void fft_rl(int, int, float *);        /* FFT function prototype */
volatile int *INTVEC = (volatile int *) 0x00; /*Address of interrupt vectors */

void recording(int *source, int rec_cmd)  /* Recording function provides I/O */
{                                         /* between the EVM and the HOST    */
  asm("  AND    0FFF8h,IF");
  asm("  OR     @_dma_int2,IE");
  dma[SOURCE] = (int)source;
  dma[DEST] = (int)host;
  dma[TRANSFER] = buffer;
  dma[GLOBAL] = 0x0E13;
  *host = rec_cmd;
}

void c_int11()    /* C_INT11 Interrupt handler function */
{
  while(*host & 0xFF);
  *host = NONE;
  dma[GLOBAL] = 0;
  asm("  AND    0FFFFh,IE");
}

void c_int99()   /* C_INT99 Interrupt handler function */
{
  for (;;);
}

void c_int05()   /* C_INT05 Interrupt handler function */
{
 PBASE[0x48] = ((int)(IO_buffer[index])) << 2;          /* Output data */
 IO_buffer[index] = (float)(PBASE[0x4C] << 16 >> 18);   /* Input data */
 if (++index >= N) index = 0;         /* Increment index, reset = N */
}

main()
{
 int loop, xfer, i;                         /* Declare variables */
 float real, img;
 INTVEC[5] = (volatile int) c_int05;     /*Install interrupt function */

 AICSET_I();                                /* Configure AIC for interrupt */
 IO_buffer = (float *) calloc(N, sizeof(float));   /* Input_Output buffer */
 data = (float *) calloc(N, sizeof(float));        /* FFT data buffer */
 while(1)
 {
   while(dma[TRANSFER]);
   do
   {
        fft_rl(N, M, (float *)data);             /* Call FFT function */
        data[0] = sqrt(data[0]*data[0])/N;       /* Magnitude of X(0) */
        for (loop = 1; loop < N/2; loop++)       /* Calculate X(1)..X(N/2-1) */
                {
                real = data[loop];               /* Real part */
                img = data[N-loop];              /* Imaginary part */
                data[loop] = sqrt(real*real+img*img)/N;    /* Find magnitude */
                }
        data[N/2] = sqrt(data[N/2]*data[N/2])/N; /* Magnitude of X(N/2) */
        for (loop = N/2+1; loop < N; loop++)     /* X(N/2+1).. X(N-1) */
                data[loop] = data[N-loop];       /* Use principle of symmetry */
        while (index);              /* Wait until the IO_buffer is empty */
        temp = data;                /* Temp pointer => Data buffer */
        data = IO_buffer;           /* IO_buffer becomes new data buffer */
        IO_buffer = temp;           /* Data buffer becomes new IO_buffer */
        IO_buffer[0] = -2048;       /* Frame sync pulse (negative spike) */
        for (xfer=0; xfer<N; xfer++)      /* Build array of N FFT values */
                fft_val[xfer] = (int)IO_buffer[xfer];

        asm("  LDI  2000h,ST");
        *host = 0x00;
        i = *host & 0xFF;
        if (i == SEND)   /* Send all N FFT values to the host for processing */
           recording(fft_val, SEND);
   }
   while (i != 128);
 }
}


