/*SUBHOST.C HOST PROGRAM FOR TRACKING OF SUBMERSIBLE*/
#include <stdlib.h>      /* Standard header files called for BC++ */
#include <stdio.h>
#include <conio.h>
#include <graphics.h>
#include "pc.c"          /* Files for EVM initialization and I/O functions */
#include "pc_1.h"        /* EVM Module PC Host Support File */
#define SEND 128
int index1 = 0;
int buffer = 512;        /* Number of FFT values sent over at one time */
int max_mag;             /* Stores Maximum Magnitude of the FFT */
int new_count;
int transfer(count)      /* Function which receives the FFT values */
int count;               /* Stores FFT point for the maximum magnitude */
{
  int i;
  int temp[512];         /* Array stores N points of the FFT */
  max_mag = 0;
  count = 0;             /* Holds FFT value where maximum magnitude occurs */
  while(READ_CMD);       /* Data is read in from EVM, using DMA transfer */
  WRITE_CMD(128);
  while (READ_CMD != 128);
  CLR_READ_ACK;
  READ_DATA;
  for(i = 0; i < buffer; i++)    /* Read takes place in a 512 byte chunk */
  {
    do
      UPDATE_STATUS0;
    while(!IS_READ_ACK);
    CLR_READ_ACK;
    temp[i] = READ_DATA;
	if (temp[i] > max_mag)       /* Maximum Magnitude is determined */
	{
	max_mag = temp[i];
	count = i;                   /* FFT value for max_mag is saved */
	}
  }
  WRITE_CMD(NONE);
  return(count);                 /* FFT value for max_mag is passed to main */
}

void main()
{
  char string[10] = {" "};       /* Graphics initialization */
  int XLOC,YLOC,XLOC1,YLOC1;
  int m,loop,loop1,fft_freq,HEAD;
  int X0 = 45, Y0 = 50;
  int X1 = X0 + 450, Y1 = Y0 + 400;
  int X10 = 550, Y10 = 50;
  int X11 = X10 + 50, Y11 = Y10 + 400;
  int maxx, maxy;
  int gdriver = DETECT, gmode, errorcode;  /* Error trapping for graphics */
  clrscr();                                /* Requires EGAVGA.BGI driver */
  HEAD = 0;
  initgraph(&gdriver, &gmode, "");
  errorcode = graphresult();
  if (errorcode != grOk)
  {
     printf("Graphics error: %s\n", grapherrormsg(errorcode));
     printf("Press any key to halt:");
     getch();
     exit(1);
   }
   setfillstyle(SOLID_FILL,BLUE);           /* Background color is selected */
   setviewport(0, 0, 700, 600, 1);
   floodfill(0, 0, BLUE);
   setbkcolor(CYAN);                        /* Geographic Display is drawn  */
   setcolor(WHITE);                         /* along with the depth gauge;  */
   for (loop = 0; loop < 21; loop++)        /* titles are entered.          */
   {
     int vert;
     vert = Y10+loop*20;
     itoa(100+(loop*20), string, 10);
     line (X10-5, vert, X10, vert);
     outtextxy (X10-30, vert, string);
   }
   settextjustify (CENTER_TEXT, CENTER_TEXT);
   settextstyle(DEFAULT_FONT, HORIZ_DIR, 2);
   outtextxy ((getmaxx() / 3)-20, Y0-30, "ANGULAR DISPLAY");
   outtextxy (X10+28, Y0-30, "DEPTH");
   setfillstyle(SOLID_FILL,CYAN);
   rectangle(370, 5, 450, 35);
   rectangle(X0-1, Y0-1, X1+1, Y1+1);
   rectangle(X10-1, Y10-1, X11+1, Y11+1);
   setviewport(X10, Y10, X11, Y11, 1);
   floodfill(0, 0, 3);
   setviewport(X0, Y0, X1, Y1, 1);
   floodfill(X0, Y0-50, 3);
   setviewport(371, 6, 449, 34, 1);
   floodfill(0, 0, 3);
   itoa(HEAD, string, 10);
   outtextxy(40, 15, string);
   XLOC = X0 +180;
   YLOC = Y0 +340;
   XLOC1 = 20;
   YLOC1 = 0;
   setviewport(X10, Y10, X11, Y11, 1);
   outtextxy (XLOC1, 105, "*");
   while (1)
	{
	int y_pos;

	delay (700);   /* 700ms delay is set for easy viewing */
	fft_freq = transfer(new_count); /* Transfer function is called and the */
	y_pos = YLOC1+105;               /* FFT value for the maximum magnitude */
	setcolor(WHITE);                /* is stored in fft_freq               */
	setviewport(X0, Y0, X1, Y1, 1);
	outtextxy (XLOC, YLOC, ".");

	if ((fft_freq > 185) && (fft_freq < 225) && (YLOC1 < 280)) /* ~4K */
		{
		setviewport(X10, Y10, X11, Y11, 1); /* If the frequency */
		setcolor(CYAN);                     /* is about 4K then the */
		outtextxy (XLOC1, y_pos, "*");   /* torpedo "dives" */
		setcolor(WHITE);                    /* shown by depth gauge */
		outtextxy (XLOC1, y_pos+20, "*");
		YLOC1 = YLOC1 + 20;
		}

	if ((fft_freq > 134) && (fft_freq < 174) && (YLOC1+80 > 0)) /* ~3K */
		{
		setviewport(X10, Y10, X11, Y11, 1); /* If the frequency */
		setcolor(CYAN);                     /* is about 3K then the */
		outtextxy (XLOC1, y_pos, "*");   /* torpedo "rises" */
		setcolor(WHITE);                    /* shown by depth gauge */
		outtextxy (XLOC1, y_pos-20, "*");
		YLOC1 = YLOC1 - 20;
		}

	if ((fft_freq > 31) && (fft_freq < 71))              /* ~1K */
		{                           /* If the frequency is */
		XLOC = XLOC + 10;           /* about 1K then the */
		HEAD = HEAD + 5;            /* torpedo steers right */
		if (HEAD > 360)
			HEAD = HEAD-360;    /* Determine and Display */
		itoa(HEAD,string,10);       /* heading of torpedo */
		setviewport(371, 6, 449, 34, 1);
		clearviewport();
		outtextxy(40, 15, string);
		}
	if ((fft_freq > 82) && (fft_freq < 122))             /* ~2K */
		{                           /* If the frequency is */
		XLOC = XLOC - 10;           /* about 2K then the */
		HEAD = HEAD - 5;            /* torpedo steers left */
		if (HEAD < 0)
			HEAD = 350-HEAD;    /* Determine and Display */
		itoa(HEAD,string,10);       /* heading of torpedo */
		setviewport(371, 6, 449, 34, 1);
		clearviewport();
		outtextxy(40, 15, string);
		}

	YLOC = YLOC -10;
	if ((XLOC >= X1-50) || (XLOC <= X0-50) || (YLOC <= Y0-60))
		{
		setviewport(X0, Y0, X1, Y1,1);
		clearviewport();
		XLOC = X0+180 ;
		YLOC = Y0 +340;
		}

	if (kbhit())                /* Any keystroke exits the program */
		{
		closegraph();
		exit(EXIT_SUCCESS);
		break;
		}
	}
}


