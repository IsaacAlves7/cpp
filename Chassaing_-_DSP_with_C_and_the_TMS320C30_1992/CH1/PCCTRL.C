/*PCCTRL.C-PC HOST INTERFACING FOR LOOP PROGRAM*/
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include "loopcmd.h"
#include "text.c"
#define READ_CMD        inp(COM_CMD)
#define WRITE_CMD(x)    outp(COM_CMD, x)
#define iobase          0x240
#define COM_CMD         iobase + 0x0800

void send_command(int i)
{
  while (READ_CMD);
  WRITE_CMD(i);
  while (READ_CMD != i);
  WRITE_CMD(0x00);
}

int get_command()
{
  unsigned int key = bioskey(0);
  fflush (stdin);
  if (key & 0x0FF) return(0x00);
  else return (key >> 8);
}

void command_process()
{
  int command;
  command = get_command();
  switch(command)
  {
    case HIGHER:
    case LOWER:
      break;
    case QUIT:
      _setcursortype(_NORMALCURSOR);
      textcolor(WHITE);
      textbackground(BLACK);
      window(1,1,80,25);
      clrscr();
      exit(EXIT_SUCCESS);
      break;
    default:
      command = 0x00;
      break;
  }
  if (command) send_command(command);
}

void menu()
{
  clrscr();
  gotoxy (1,2);
  printf("  F5    DECREASE GAIN");
  gotoxy (1,3);
  printf("  F6    INCREASE GAIN");
  gotoxy (1,4);
  printf("  END   QUIT");
}

void graphics()
{
  _setcursortype(_NOCURSOR);
  textcolor (LIGHTGRAY);
  textbackground (BLUE);
  clrscr();
  shadow_box (25,9,55,15,LIGHTGRAY, BLUE);
  window (26,10,54,14);
}

void main()
{
  clrscr();
  graphics();
  for (;;)
  {
    menu();
    command_process();
  }
}





