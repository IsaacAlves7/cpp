/*VIDEOPC.C-PC HOST PROGRAM FOR IMAGE PROCESSING PROJECT*/
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <graphics.h>
#include "pc.c"
#include "pc_1.h"
#include "projcmd.h"
void send_cmd();
void init_graphic();
void transfer();
int pause = 0;
int update = 1;

void main()
{
  init_evm();
  init_graphic();
  buffer = 256;
  send_cmd();
  for(;;)
  {
    if (!pause) transfer();
    if (kbhit()) update = 1;
    if (kbhit()) send_cmd();
  }
}

void send_cmd()
{
  int command;
  command = get_command();
  switch(command)
  {
    case RESET:
    case LPA:
    case LPB:
    case LPC:
    case SEL_RATE:
    case INCREASE:
    case DECREASE:
    case AVER:
    case EDGE1:
    case EDGE2:
    case L_PLUS:
    case L_MINUS:
	 break;
    case PAUSE: if (pause) {pause = 0; break;}
		else {pause = 1; break;}
    case QUIT:
	 closegraph();
	 exit(EXIT_SUCCESS);
	 break;
    default:
	 command = NONE;
	 break;
  }
  if (command) send_command(command);
}

void init_graphic()
{
  int m, loop;
  int X0 = 65, Y0 = 50;
  int X1 = X0 + 508, Y1 = Y0 + 265;
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
   for (loop = 0; loop < 11; loop++)
   {
     int vert;
     char string[4] = {0};
     vert = Y0+5+loop*25;
     itoa(100-(loop*10), string, 10);
     line (X0-5, vert, X0, vert);
     outtextxy (X0-20, vert, string);
   }
   settextstyle(DEFAULT_FONT, HORIZ_DIR, 2);
   outtextxy (getmaxx()/2, Y0-30, "VIDEO LINE RATE PROCESSING");
   settextstyle(DEFAULT_FONT, VERT_DIR, 1);
   outtextxy (X0-50, Y0+128, "IRE units");
   settextstyle(DEFAULT_FONT, HORIZ_DIR, 0);
   settextjustify(LEFT_TEXT, CENTER_TEXT);
   outtextxy (94, 360, " F1  RESET    ");
   outtextxy (94, 375, " F2  500K LPF ");
   outtextxy (94, 390, " F3  1Meg LPF ");
   outtextxy (94, 405, " F4  3Meg LPF ");
   outtextxy (94, 420, " F5  AVERAGING");
   outtextxy (94, 442, "     LINE SEL ");
   outtextxy (360,360, " F6  # LINES  ");
   outtextxy (360,375, " F7  EDGE1    ");
   outtextxy (360,390, " F8  EDGE2    ");
   outtextxy (360,405, " F9  PAUSE    ");
   outtextxy (360,420, "F10  QUIT     ");
   outtextxy (360,442, "     RATE     ");
   line (108,432,111,432);                              // start of up arrow
   line (107,433,112,433);
   line (106,434,113,434);
   line (109,431,109,439);
   line (110,431,110,439);
   line (109,446,109,454);                              // start of down arrow
   line (110,446,110,454);
   line (106,451,113,451);
   line (107,452,112,452);
   line (108,453,111,453);
   line (372,435,380,435);                            // start of left arrow
   line (372,436,380,436);
   line (375,439,375,432);
   line (374,438,374,433);
   line (373,437,373,434);
   line (372,450,380,450);                            // start of right arrow
   line (372,451,380,451);
   line (377,454,377,447);
   line (378,453,378,448);
   line (379,452,379,449);
   rectangle (210,350,285,460);
   rectangle (470,350,545,460);
   setfillstyle(SOLID_FILL,LIGHTGRAY);
   setviewport(X0, Y0, X1, Y1, 1);
   outtextxy(150,128, "PRESS SPACE BAR TO CONTINUE");
}

void transfer()
{
  int i, value, index, status, line, sel_rate, rate;
  char string[4] = {0};
  int temp[720];
  while(READ_CMD);
  WRITE_CMD(128);
  while(READ_CMD != 128);
  CLR_READ_ACK;
  READ_DATA;
  index = 0;
  for(i = 0; i < buffer; i++)
  {
    do
    {
      UPDATE_STATUS0;
    }
    while(!IS_READ_ACK);
    CLR_READ_ACK;
    value = (int)READ_DATA;
    temp[index++] = (value >> 8) & 255;
    temp[index++] = value & 255;
  }
  WRITE_CMD(NONE);
  bar(0,0,508,265);
  setcolor(LIGHTGRAY);
  lineto (0, (256- temp[0]));
  setcolor(WHITE);
  for (i = 0; i < 508; i++)
    lineto (i, (256 - temp[i]));
  if (update)
  {
    update = 0;
    status = temp[508];
    line = temp[509];
    rate = temp[510];
    sel_rate = temp[511];
    setviewport(211,351,284,459,1);
    bar(0,0,73,108);
    if (status & 1) outtextxy(25, 24,"ON");
    else outtextxy(25,24, "OFF");
    status = status >> 1;
    if (status & 1) outtextxy(25, 39,"ON");
    else outtextxy(25,39, "OFF");
    status = status >> 1;
    if (status & 1) outtextxy(25, 54,"ON");
    else outtextxy(25,54, "OFF");
    status = status >> 1;
    if (status & 1) outtextxy(25, 69,"ON");
    else outtextxy(25,69, "OFF");
    itoa(line, string, 10);
    outtextxy (25,91, string);
    setviewport(471,351,544,459,1);
    bar (0,0,73,108);
    itoa(sel_rate, string, 10);
    outtextxy (33, 9, string);
    status = status >> 1;
    if (status & 1) outtextxy(25, 24,"ON");
    else outtextxy(25, 24, "OFF");
    status = status >> 1;
    if (status & 1) outtextxy(25, 39,"ON");
    else outtextxy(25, 39, "OFF");
    rate = 60/rate;
    if (rate > 1)
    {
      itoa(rate, string, 10);
      outtextxy (25, 92, "1/  ");
      outtextxy (45, 92, string);
    }
    else
    {
      rate = temp[510];
      rate = rate/60;
      if (!rate) rate = 1;
      itoa(rate, string, 10);
      outtextxy (33, 92, string);
    }
    setviewport(65,50,573,315,1);
  }
}







