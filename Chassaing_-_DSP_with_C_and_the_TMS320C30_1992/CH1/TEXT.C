/*TEXT.C- C functions for text displays */
#include <dos.h>

void writechar (char ch, int back_color, int fore_color, int repeat_char);

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

void writechar (char ch, int back_color, int fore_color, int repeat_char)
{
   int att;
   back_color <<= 4;                      //shift color attr backgnd bits
   att = (back_color + fore_color);       // add backgnd and foregnd attr bits
   _AH = 9;                               // interrupt 0x10 sub-function 9
   _AL = ch;                              // character to be output 
   _BH = 0;                               // video page 
   _BL = att;                             // video attribute 
   _CX = repeat_char;                     // repetition factor 
   geninterrupt(0x10);                    // output the char 
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
