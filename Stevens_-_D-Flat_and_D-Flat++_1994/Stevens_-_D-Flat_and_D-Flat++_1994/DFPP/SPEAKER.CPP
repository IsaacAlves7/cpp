// -------- speaker.cpp

#include <dos.h>
#include <conio.h>
#include "speaker.h"

// -------- sound a tone
void Speaker::Beep()
{
    outp(0x43, 0xb6);               // program the frequency
    outp(0x42, (int) (COUNT % 256));
    outp(0x42, (int) (COUNT / 256));
    outp(0x61, inp(0x61) | 3);      // start the sound
    // -------- wait two clock ticks
    const int far *clk = (int far *) MK_FP(0x40,0x6c);
	int then = *clk+2;
	while (*clk < then)
		;
    outp(0x61, inp(0x61) & ~3);     // stop the sound
}

