/*PCLOOP.C COMMAND PROGRAM FOR COMMUNICATION BETWEEN PC/EVM */
#include <stdio.h>
#include <dos.h>

#define IOBASE          0x240                 // base address on PC
#define COM_CMD         IOBASE + 0x0800       // addr of 8 bit R/W port

void main()
{
  int i;                                      // declare variable
  clrscr();                                   // clear screen
  outp(COM_CMD,0x00);                         // clear 8 bit R/W port to zero
  for (;;)                                    // create endless loop
  {
    do
    {
      printf("\nPRESS CTRL BREAK TO QUIT");
      printf("\nEnter attenuation value (1-10) : ");
      scanf("%i", &i);                       // get attenuation value
    }
    while ((i<1) || (i>10));                 // repeat if value not valid
    outp(COM_CMD,i);                         // output value to R/W port
    clrscr();                                // clear screen for next value
  }
}





