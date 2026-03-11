/*DMAHOST.C-HOST COMMUNICATION WITH DMA*/
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include "pc.c"
#include "pc_1.h"
#define SEND 128
void transfer();
int index1 = 0;
int data[10] = {1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};
int buffer = 10;

void main()
{
  int i;
  init_evm();
  while(READ_CMD);
  WRITE_CMD(64);
  while(READ_CMD != 64);
  CLR_WRITE_ACK;
  for (i = 0; i < buffer; i++)
  {
    WRITE_DATA(data[i]);
    do
      UPDATE_STATUS0;
    while (!IS_WRITE_ACK);
    CLR_WRITE_ACK;
  }
  WRITE_CMD(NONE);
  transfer();
}

void transfer()
{
  int i;
  int temp[10];

  while(READ_CMD);
  WRITE_CMD(128);
  while (READ_CMD != 128);
  CLR_READ_ACK;
  READ_DATA;
  for(i = 0; i < buffer; i++)
  {
    do
      UPDATE_STATUS0;
    while(!IS_READ_ACK);
    CLR_READ_ACK;
    temp[i] = READ_DATA;
    printf("%i\n", temp[i]);
  }
  WRITE_CMD(NONE);
}







