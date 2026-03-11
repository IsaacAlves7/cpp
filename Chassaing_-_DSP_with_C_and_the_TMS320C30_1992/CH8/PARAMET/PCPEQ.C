/*PCPEQ.C-PC HOST/EVM COMMUNICATION PROGRAM*/
#include <stdio.h>
#include <dos.h>
#define IOBASE 0x280
#define PORT (IOBASE+0x800)
main()
    {
    int i;
    outp(PORT,0x00);
    clrscr();
    do
	{
	printf("\nEnter filter number (1-16) : ");
	scanf("%i",&i);
	}
    while ((i<1) || (i>16));
    outp(PORT,i);
    }

