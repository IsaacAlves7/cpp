/*AMPLIT.CPP-TO PLOT AMPLITUDE AND PHASE
This program was compiled with Turbo C++. Use the large memory model specified 
in the compiler option. This program uses interrupts to read and write to the
display by calling bios routines. After the data is entered the information is
read out of display memory and placed into a string.  The string is converted
to a float where the information is used to determine the transfer function of
the filter. The magnitude and the phase of the transfer function is plotted to
the graphics display.By pressing F1 the graph is read from the display and sent
to a laser printer where a hard copy of the screen is obtained. */
#include <math.h>
#include <complex.h>
#include <dos.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <graphics.h>
#include <string.h>
#define pi 3.1415926
void writechar (char ch, int back_color, int fore_color, int repeat_char);
/*draws a box at defined coordinates and adds a black shadow to the bottom 
  and right of the box */
void shadow_box (X1, Y1, X2, Y2, back_color, fore_color)
{
  int loop;
  if (X1 < X2 && X2 <= 80 && Y1 < Y2 && Y2 <= 25)
    {
      window (X1+2,Y1+1,X2+2,Y2+1);
      textattr ((BLACK << 4) + BLACK);
      clrscr();
      window (X1,Y1,X2,Y2);
      textattr ((back_color << 4) + fore_color);
      clrscr();
      window (1,1,80,25);
      gotoxy (X1,Y1);
      writechar (201, back_color, fore_color, 1);
      gotoxy (X1 + 1, Y1);
      writechar (205, back_color, fore_color, (X2 - X1) - 1 );
      gotoxy (X2,Y1);
      writechar (187, back_color, fore_color, 1);
      for (loop = Y2 - 1; loop > Y1; --loop)
	{
	  gotoxy (X1,loop);
	  writechar (186, back_color, fore_color, 1);
	  gotoxy (X2,loop);
	  writechar (186, back_color, fore_color, 1);
	}
      gotoxy (X1,Y2);
      writechar (200, back_color, fore_color, 1);
      gotoxy (X1+1,Y2);
      writechar (205, back_color, fore_color, (X2 - X1) - 1 );
      gotoxy (X2,Y2);
      writechar (188, back_color, fore_color, 1);
      gotoxy (X1+1,Y1+1);
    }
}
/* outputs a character at the current cursor position using the video BIOS
   to avoid the scrolling of the screen when writing to location (80,25).  */

void writechar (char ch, int back_color, int fore_color, int repeat_char)
{  int att;
   back_color <<= 4;                      //shift color attr backgnd bits
   att = (back_color + fore_color);       // add backgnd and foregnd attr bits
   _AH = 9;                               // interrupt 0x10 sub-function 9
   _AL = ch;                              // character to be output 
   _BH = 0;                               // video page 
   _BL = att;	                          // video attribute 
   _CX = repeat_char;                     // repetition factor 
   geninterrupt(0x10);                    // output the char 
}

void LJ_Graphic()
{
  int xaspect, yaspect, maxX, maxY, line, xword, pixel, xwidth, ywidth;
  double xpos, ypos, prnstep, ratio;
  char chr;
				   
  maxX = getmaxx();                         //get number of horizontal pixel
  maxY = getmaxy();                         //get number of vertical pixel
  getaspectratio(&xaspect, &yaspect);       //get the screen aspect ratio
  ratio = (double) xaspect/ (double) yaspect;      
  setviewport(0,0,maxX,maxY,0);             //set viewport for full screen
  xpos = 690;                               //initial position of prn cursor
  ypos = 500;                               //initial position of prn cursor
  prnstep = 7.2/ratio;                      //match printer aspr fo screen aspr
  fprintf (stdprn, "\x1B&E\x1B&11H\x1B&1O\x1B*p0X\x1B*p0Y\x1B*t100R");
  for (line = 0; line <= maxY; line++)
  {
    ywidth = 6;
    if (ypos < 1000.0) ywidth--;
    if (ypos < 100.0) ywidth--;
    if (ypos < 10.0) ywidth--;
    fprintf (stdprn, "\x1B&a%-*.1fh%-*.1fV", 5, xpos, ywidth, ypos);
    ypos += prnstep;
    fprintf (stdprn, "\x1B*r1A\x1B*b%dW", maxX/8);
    for (xword = 0; xword < maxX/8; xword++)
    {
      chr = 0;
      for (pixel = 0; pixel < 8; pixel++)  
      {                     //reads series of 8 pixels to create graphics char
	chr <<= 1;
	if (getpixel (xword*8+pixel, line)) chr++;
      }
      fprintf (stdprn, "%c", chr);             //sends graphic char to printer
    }
    fprintf (stdprn, "\x1B*rB");               //ends graphics line
  }
  fprintf (stdprn, "\x0C\x1B&10O\x1B&11H\x1B&E");
}

/* reads the character from display memory at current cursor position */
char read_char (void)     
{
  char ch;
  _AH = 8;                      // call video service 8
  _BH = 0;                      // set display page to 0
  geninterrupt(0x10);           // call video interrupt
  ch = _AL;                     // screen char in AL
  return(ch);
}

void reverse_video (void)
{
  textcolor(BLUE);
  textbackground(LIGHTGRAY);
}

void normal_video (void)
{
  textcolor(LIGHTGRAY);
  textbackground(BLUE);
}

/* prepares the display to switch the highlighted line to new position */
void switch_line (void)
{
   int loop;
   char str[10];
   gotoxy(1,wherey());
   for (loop = 0; loop <= 9; loop++)
   {
     str[loop] = read_char();
     gotoxy(wherex()+1,wherey());
   }
   normal_video();
   gotoxy(1,wherey());
   for (loop = 0; loop <= 9; loop++)
     putch(str[loop]);
}

void help (void)
{
  struct text_info textinfo;
  int top = 4;
  int bottom = 19;
  int left = 5;
  int right = 75;
  char buffer[4000];
  char ch;
  gettextinfo(&textinfo);
  window(1,1,80,25);
  gettext(left,top,right+2,bottom+1,buffer);
  shadow_box(left,top,right,bottom,CYAN,WHITE);
  window(left+1,top+1,right-1,bottom-1);
  clrscr();
  gotoxy (30,1);
  cprintf ("HELP SCREEN");
  gotoxy (2,3);
  textcolor(BLACK);
  cprintf ("Input the coefficients to the appropriate powers of z");
  gotoxy (2,4);
  cprintf ("The following keys can be used to enter the coefficients");
  textcolor(WHITE);
  gotoxy (10,6);
  cprintf ("0 - 9");
  gotoxy (10,7);
  cprintf ("+ - .");
  gotoxy (5,8);
  cprintf ("ARROW KEYS");
  gotoxy (8,9);
  cprintf ("TAB KEY");
  gotoxy (6,10);
  cprintf ("ENTER KEY");
  gotoxy (5,11);
  cprintf ("DELETE KEY");
  gotoxy (6,12);
  cprintf ("BACKSPACE");
  gotoxy (23,14);
  _setcursortype(_NOCURSOR);
  cprintf ("Press any KEY to continue");
  textcolor(BLACK);
  gotoxy(17,6);
  cprintf ("Enters numeric values in the field");
  gotoxy(17,7);
  cprintf ("Enters the sign or decimal point of the value");
  gotoxy(17,8);
  cprintf ("Moves the cursor in the corresponding direction");
  gotoxy(17,9);
  cprintf ("Moves cursor from numerator to denominator");
  gotoxy(17,10);
  cprintf ("Moves cursor to the begining of the next field");
  gotoxy(17,11);
  cprintf ("Deletes the entry at the current cursor position");
  gotoxy(17,12);
  cprintf ("Moves cursor one position left and deletes entry");
  while(!kbhit());
  ch = getch();
  if (ch == 0) 
    ch = getch();
  _setcursortype(_NORMALCURSOR);
  puttext(left,top,right+2,bottom+1,buffer);
  window(textinfo.winleft,textinfo.wintop,textinfo.winright,textinfo.winbottom);
  gotoxy(textinfo.curx,textinfo.cury);
  textattr(textinfo.attribute);
}
void beep (void)
{
  sound(500);
  delay(100);
  nosound();
}
	
/* This function creates a data entry box and several editing keystrokes
   to enable the user to easily input his information in the proper area.
   The function also limits the available keystrokes to those required
   for data entry.  The user moves the cursor around the data entry box
   and is finished when F10 is pressed.  */ 

void coeff_box (void)
  {
    int test, loop, tab, num = 1, X, period;
    char str[10];
    char ch, fk = 0;
    textbackground(LIGHTGRAY);
    clrscr();
    shadow_box (20,3,60,21,BLUE,WHITE);
    window(21,4,59,20);
    clrscr();
    window(20,3,60,21);
    gotoxy(11,2);
    cprintf ("FILTER COEFFICIENTS");
    gotoxy(6,18);
    cprintf ("F1 HELP    F5 Quit    F10 PLOT");
    gotoxy(7,4);
    cprintf ("NUMERATOR");
    gotoxy(26,4);
    cprintf ("DENOMINATOR");
    gotoxy(1,3);
    writechar(204,BLUE,WHITE,1);
    gotoxy(2,3);
    writechar(205,BLUE,WHITE,39);
    gotoxy(41,3);
    writechar(185,BLUE,WHITE,1);
    gotoxy(1,5);
    writechar(199,BLUE,WHITE,1);
    gotoxy(1,17);
    writechar(199,BLUE,WHITE,1);
    gotoxy(2,5);
    writechar(196,BLUE,WHITE,39);
    gotoxy(41,5);
    writechar(182,BLUE,WHITE,1);
    gotoxy(41,17);
    writechar(182,BLUE,WHITE,1);
    gotoxy(2,17);
    writechar(196,BLUE,WHITE,39);
    gotoxy(20,3);
    writechar(209,BLUE,WHITE,1);
    gotoxy(20,17);
    writechar(193,BLUE,WHITE,1);
    for (loop = 4; loop <= 16; loop++)
    {  
      gotoxy(20,loop);
      writechar (179,BLUE,WHITE,1);
    }
    gotoxy(20,5);
    writechar (197,BLUE,WHITE,1);
    for (loop = 6; loop <= 15; loop++)
    {  
      gotoxy(3,loop);
      cprintf("z-%i", loop-6);
      gotoxy(22,loop);
      cprintf("z-%i", loop-6);
    }
    gotoxy(2,16);
    cprintf("z-10");
    gotoxy(21,16);
    cprintf("z-10");
    window(27,8,59,18);
    gotoxy(1,1); 
    textcolor (LIGHTGRAY);
    test = 0;
    do                                                 // data entry loop
    {
      tab = 0;
      for (loop = 0; loop <= 9; loop++)
      {
	str[loop] = read_char();
        gotoxy(wherex()+1,wherey());
      }
      reverse_video();
      gotoxy(1,wherey());
      for (loop = 0; loop <= 9; loop++)               //creates reversed bar
      { 
	putch(str[loop]);
      }
      gotoxy(1,wherey());
      do
	{
	  tab = 0;
	  ch = getch();
	  if (ch >= '0' && ch <= '9')
	  {
	    if (wherex() != 10)
	      putch(ch);
	    else
	    {   
	      writechar(ch,LIGHTGRAY,BLUE,1);
	      beep();
	    }
	  }
	  else 
	    if (ch != 0)
	    { 
	      switch(ch)
	      {
		case 13 : fk = 1; ch = 80; break;    //enter key pressed
		case 8 :                             //backspace key
		  if (wherex() != 1)
		  {
		    writechar(0,LIGHTGRAY,BLUE,1);
		    putch(ch);
		    writechar(0,LIGHTGRAY,BLUE,1);
                  }
		  else beep();
		  break;
		case '.' :
		  X = wherex();
		  gotoxy(1, wherey());
		  period = 0;
		  for (loop = 1; loop <= 10; loop++)
		  {
		    gotoxy(loop, wherey());
		    if (read_char() == '.') period = 1;
                  }
		  gotoxy(X,wherey());
		  if (period == 1) beep();
		  else
		    putch(ch); 
		  break;
		case '+' : 
		  if (wherex() == 1) 
		    putch(ch); 
		  else
		    beep();
		   break;
		case '-' : 
		  if (wherex() == 1) 
		    putch(ch); 
		  else
		    beep();
		  break;
		case 9 : 
		  if (num)
		    num = 0;
		  else
		    num = 1;
		  switch_line();
		  if (!num)
		    window(47,8,59,18);
		  else
		    window(27,8,59,18);
		  gotoxy(1,1);
                  for (loop = 0; loop <= 9; loop++)
                  {
	            str[loop] = read_char();
                    gotoxy(wherex()+1,wherey());
                  }
                  reverse_video();
                  gotoxy(1,wherey());
                  for (loop = 0; loop <= 9; loop++)
                  { 
	            putch(str[loop]);
                  }
                  gotoxy(1,wherey());
		  tab = 1;
                break;
		default : beep(); break;
	      }
            }
	    else
	    {
	      ch = getch();
              fk = 0;
	      switch (ch)                          
	      {
		 case 75 :                         //left arrow key pressed
		   if (wherex() != 1) 
		     gotoxy(wherex()-1,wherey());
		   else
		     beep();
		   break;
		 case 77 :                         //right arrow key pressed
		   if (wherex() != 10)
		     gotoxy(wherex()+1,wherey());
		   else
		     beep();
		   break;
		 case 59 : help(); break;          // F1 key pressed
		 case 63 :                         // F5 key pressed
		   window (1,1,80,25);
		   textbackground(BLACK);
		   textcolor(WHITE);
		   clrscr(); 
		   exit(0);                 
		 default : fk = 1;
	       }
	     }
      }
      while (!fk || tab);
      if (fk == 1)
      {
	 fk = 0;
	 switch (ch)
	 {
           case 72 :                          //up arrow key pressed
	     switch_line();
	     if (wherey() != 1) 
	       gotoxy(1,wherey()-1);
	     else
	       gotoxy(1,11);
	     break;
	   case 80 :                         //down arrow key pressed
	     switch_line();
	     if (wherey() != 11)
	       gotoxy(1,wherey()+1);
	     else
	       gotoxy(1,1);
	     break;
	   case 68 : test = 1; break;       // F10 key pressed
	   default : beep(); gotoxy(1,wherey()); break;
	 }
      }
    }
    while (!test);
  }

main()
{
  int left, top, right, bottom, maxx, maxy, loop, loop1;
  double v, mag, mag_1, real_part, imag_part, phase, phase_1;
  char str[10];
  char ch;
  double A[11] = {0.0};
  double B[11] = {0.0};
  complex NUM, DEN;
     textmode(C80);  
  coeff_box();
  window(27,8,59,18);
  for (loop = 0; loop <= 10; loop++)        //reads data info from the display
  {
    for (loop1 = 1; loop1 <= 10; loop1++)
    {
      gotoxy(loop1,loop+1);                 //move cursor to read next char
      str[loop1-1] = read_char();
    } 
    A[loop] = atof(str);                    //converts string to float
    for (loop1 = 1; loop1 <= 10; loop1++)
    {
      gotoxy(loop1+20,loop+1);
      str[loop1-1] = read_char();
    }
    B[loop] = atof(str);
  }
  window(1,1,80,25);
  clrscr();
  
  int gdriver = DETECT, gmode, errorcode;    // request auto detection
  initgraph(&gdriver, &gmode, "");           // initialize graphics mode
  errorcode = graphresult();                 // read result of initialization 
  if (errorcode != grOk)                     // an error occurred
    {
       printf("Graphics error: %s\n", grapherrormsg(errorcode));
       printf("Press any key to halt:");
       getch();
    }

  maxx = getmaxx();                       //number of horizonal screen pixel
  maxy = getmaxy();                       //number of vertical screen pixel
  left = maxx/2 - 200;                    //left side of graph
  right = maxx/2 + 200;                   //right side of graph
  top = maxy/2 - 170;                     //top side of graph
  bottom = maxy/2 + 130;                  //bottom side of graph
  rectangle (left, top, right, bottom);
  setviewport (left, top, right, bottom, 1);
  for (loop = 1; loop <= 5; loop++)
  {
    line(loop*80,300,loop*80,292);        //draws vert hash marks
    line(0,loop*100,5,loop*100);          //draws horz hash marks
  }
  setlinestyle (SOLID_LINE,1,1);
  setviewport (left, top, right, bottom, 0);
  settextjustify(CENTER_TEXT,CENTER_TEXT);
  outtextxy(200,330, "Normalized Frequency = f/fN");
  outtextxy(200,345, "F1 for PRINTOUT      ENTER to continue");

  for (loop = 0; loop <= 5; loop++)
  {
    gcvt(loop/5.0, 10, str);              //converts float to string
    outtextxy(loop*80,310,str);
  }
  settextstyle(DEFAULT_FONT, VERT_DIR, 1);
  outtextxy(-75,150, "Magnitude");
  settextstyle(DEFAULT_FONT, HORIZ_DIR, 1);
  settextjustify (RIGHT_TEXT,CENTER_TEXT);
  outtextxy(-10,300,"0.0");
  outtextxy(-10,200,"0.5");
  outtextxy(-10,100,"1.0");
  outtextxy(-10,0,"1.5");
  setviewport (left, top, right, bottom, 1);
  for (loop = 1; loop <= 400; loop++)        //loop to plot graph
  {
    v = loop/400.0;
    NUM = complex(A[0],0);                   //numerator 
    DEN = complex(B[0],0);                   //denominator 
    for (loop1 = 1; loop1 <= 10; loop1++)    //transfer fuction
    {
      NUM += complex(A[loop1]*cos(pi*v*loop1),-A[loop1]*sin(pi*v*loop1));
      DEN += complex(B[loop1]*cos(pi*v*loop1),-B[loop1]*sin(pi*v*loop1));
    }
    mag = sqrt(norm(NUM/DEN));               //magnitude of transfer function
    if (loop == 1) 
      mag_1 = mag;
    else                                     //plot graph
      {
	line (loop-1, 300 - mag_1*200, loop, 300 - mag*200);
	mag_1 = mag;
      }
   } 
  while(!kbhit());
  ch = getch();
  if (ch == 0)
  {
    ch = getch();
    if (ch == 59)
      LJ_Graphic();                          //print graphics
  }
  setviewport(0,0,115,maxy,0);
  clearviewport();
  setviewport(left+1,top+1,right-1,bottom-1,1);
  clearviewport();
  setviewport(left,top,right,bottom,0);
  setlinestyle(DASHED_LINE,1,1);
  line(0,150,400,150);
  setlinestyle(SOLID_LINE,1,1);
  for (loop = 1; loop <= 5; loop++)
    line(loop*80,300,loop*80,292);        //draws vert hash marks
  outtextxy(-15,0,"1.57");
  outtextxy(-15,150,"0.0");
  outtextxy(-15,300,"-1.57");
  settextstyle(DEFAULT_FONT,VERT_DIR,1);
  outtextxy(-75,150,"RADIANS");
  settextstyle(DEFAULT_FONT,HORIZ_DIR,1);
  setviewport (left, top, right, bottom, 1);
  for (loop = 1; loop <= 400; loop++)
  {
    v = loop/400.0;
    NUM = complex(A[0],0);
    DEN = complex(B[0],0);
    for (loop1 = 1; loop1 <= 10; loop1++)
    {
      NUM += complex(A[loop1]*cos(pi*v*loop1),-A[loop1]*sin(pi*v*loop1));
      DEN += complex(B[loop1]*cos(pi*v*loop1),-B[loop1]*sin(pi*v*loop1));
    }
    real_part = real(NUM/DEN);
    imag_part = imag(NUM/DEN);
    phase = atan(imag_part/real_part);           //phase of transfer function
  if (loop == 1)
     phase_1 = phase;
  else                                           //plot phase
    {
      line (loop-1, 150-(150*phase_1/(pi*0.5)), loop, 150-(150*phase)/(pi*0.5));
      phase_1 = phase;
    }
  }          
  while (!kbhit());
  ch = getch();
  if (ch == 0)
  {
    ch = getch();
    if (ch == 59)
      LJ_Graphic();
  }
  closegraph();
  textcolor(WHITE);
  textbackground(BLACK);
  clrscr();
  return(0);
}    
  



  











