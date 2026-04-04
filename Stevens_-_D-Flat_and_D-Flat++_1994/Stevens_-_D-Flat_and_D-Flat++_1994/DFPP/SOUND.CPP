// -------- sound.cpp

#include <dos.h>
#include "sound.h"

// -------- sound a tone
void Sound::Beep()
{
    Wait(1);
    outp(0x43, 0xb6);               // program the frequency
    outp(0x42, (int) (COUNT % 256));
    outp(0x42, (int) (COUNT / 256));
    outp(0x61, inp(0x61) | 3);      // start the sound
    Wait(1);
    outp(0x61, inp(0x61) & ~3);     // stop the sound
}

// ------- wait n clock ticks
void Sound::Wait(int n)
{
	static int far *clk = (int far *) MK_FP(0x40,0x6c);
    int now = *clk;
	int then = now+n;
    while (now < then)
        now = *clk;
}


