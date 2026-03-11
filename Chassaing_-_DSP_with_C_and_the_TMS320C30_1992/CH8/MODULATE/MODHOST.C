/*MODHOST.C  PC Host/EVM Modulation Program  */
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include "pc.c"
#include "pc_1.h"
#define buffer 1          /*one DMA word*/
int c, null, input;
int dma_data[1];
main()
{
  int i, j;
  init_evm();            /* init evm for DMA Comm's */
  while(READ_CMD);       /* wait for EVM/DMA ready  */
  c = '0';
  while (c != 'x')
  {
   printf("Enter 1 to obtain a 1-kHz output\n");
   printf("Enter 2 to obtain a 2-kHz output\n");
   printf("Enter 3 to obtain a 3-kHz output\n");
   printf("Enter 4 to obtain a 4-kHz output\n\n");
   printf("Enter selection or x to exit ====>","\n");
   c = getchar();
   null = getchar();
   switch(c)
   {
	case '1':
	 dma_data[0] = 1;
	 break;
	case '2':
	 dma_data[0] = 1000;
	 break;
	case '3':
	 dma_data[0] = 2000;
	 break;
	case '4':
	 dma_data[0] = 3000;
	 break;
	default:
	 break;
   }
   WRITE_CMD(64);              /* send DMA read command */
   while(READ_CMD != 64);      /* readback command      */
   CLR_WRITE_ACK;
	 for(j=0; j < buffer; j++)
	 {
	  WRITE_DATA(dma_data[j]);   /*send EVM new frequency*/
	  do
	   UPDATE_STATUS0;
	  while (!IS_WRITE_ACK);
	  CLR_WRITE_ACK;
	  }
   WRITE_CMD(NONE);
  }
 }




