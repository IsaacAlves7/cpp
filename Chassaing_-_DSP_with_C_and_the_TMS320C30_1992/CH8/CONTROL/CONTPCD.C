/*CONTPCD.C-FOR PID CONTROLLER*/
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <graphics.h>
#include "pc.c"
#include "pc_1.h"
#define SEND 128
void init_graphic();
void transfer();
int index1 = 0;
int data[3] = {0};
float input;

void main()
{
  int i;
  init_evm();
  clrscr();
  printf("\nInput a value for PGAIN (0 - 1000)  ");
  scanf("%i", &data[0]);
  printf("\nInput a value for IGAIN (0.0001 - 1.0)  ");
  scanf("%f",  &input);
  data[1] = (int)(input*1000);
  printf("\nInput a value for DGAIN (0 - 1000)  ");
  scanf("%i", &data[2]);
  init_graphic();
  bar(0,0,512,265);
  buffer = 3;
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
  for(;;)
  {
    transfer();
    if (kbhit())
    {
      closegraph();
      exit(EXIT_SUCCESS);
      break;
    }
  }
}

void init_graphic()
{
  char string[10] = {" "};
  int m, loop;
  int X0 = 65, Y0 = 50;
  int X1 = X0 + 512, Y1 = Y0 + 265;
  int gdriver = DETECT, gmode, errorcode;
  initgraph(&gdriver, &gmode, "");
  errorcode = graphresult();
  if (errorcode != grOk)
  {
     printf("Graphics error: %s\n", grapherrormsg(errorcode));
     printf("Press any key to halt:");
     getch();
     exit(1);
   }
   setbkcolor(BLUE);
   rectangle(X0-1, Y0-1, X1+1, Y1+1);
   settextjustify (CENTER_TEXT, CENTER_TEXT);
   for (loop = 0; loop < 13; loop++)
   {
     int vert;
     vert = Y0+12+loop*20.3;
     itoa(120-(loop*10), string, 10);
     line (X0-5, vert, X0, vert);
     outtextxy (X0-20, vert, string);
   }
   settextstyle(DEFAULT_FONT, HORIZ_DIR, 2);
   outtextxy (getmaxx()/2, Y0-30, "TMS320C30 PID CONTROLLER");
   settextstyle(DEFAULT_FONT, HORIZ_DIR, 0);
   settextjustify(LEFT_TEXT, CENTER_TEXT);
   outtextxy (100, Y1+50, "PGAIN =");
   itoa(data[0], string, 10);
   outtextxy (165, Y1+50, string);
   outtextxy (100, Y1+70, "IGAIN =");
   gcvt(input, 6, string);
   outtextxy (165, Y1+70, string);
   outtextxy (100, Y1+90, "DGAIN =");
   itoa(data[2], string, 10);
   outtextxy (165, Y1+90, string);
   setfillstyle(SOLID_FILL,LIGHTGRAY);
   setviewport(X0, Y0, X1, Y1, 1);
}

void transfer()
{
  int i, index;
  int temp[6];
  char ch;
  while(READ_CMD);
  WRITE_CMD(128);
  while (READ_CMD != 128);
  CLR_READ_ACK;
  READ_DATA;
  for(i = 0; i < buffer; i++)
  {
    temp[buffer+i] = temp[i];
    do
      UPDATE_STATUS0;
    while(!IS_READ_ACK);
    CLR_READ_ACK;
    temp[i] = READ_DATA;
  }
  WRITE_CMD(NONE);
  if (index1 == 512)
  {
    bar(0,0,512,265);
    index1 = 0;
  }
  for (i = 0; i < buffer; i++)
  {
    putpixel (index1, 256 - temp[i], 3+i);
  }
  index1++;
  setviewport(65,50,577,315,1);
}





