// ----------- keyboard.cpp

#include <stdio.h>
#include <bios.h>
#include <dos.h>
#include <conio.h>
#include "desktop.h"

/* ----- table of alt keys for finding shortcut keys ----- */
static int altconvert[] = {
    ALT_A,ALT_B,ALT_C,ALT_D,ALT_E,ALT_F,ALT_G,ALT_H,
    ALT_I,ALT_J,ALT_K,ALT_L,ALT_M,ALT_N,ALT_O,ALT_P,
    ALT_Q,ALT_R,ALT_S,ALT_T,ALT_U,ALT_V,ALT_W,ALT_X,
    ALT_Y,ALT_Z,ALT_0,ALT_1,ALT_2,ALT_3,ALT_4,ALT_5,
    ALT_6,ALT_7,ALT_8,ALT_9
};

Keyboard::Keyboard()
{
	shift = GetShift();
	keydown = False;
	insert = True;
}

/* ---- Test for keystroke ---- */
Bool Keyboard::KeyHit()
{
    _AH = 1;
    geninterrupt(KEYBRD);
    return (Bool)((_FLAGS & ZEROFLAG) == 0);
}

/* ---- Read a keystroke ---- */
int Keyboard::GetKey()
{
    int c;
    while (KeyHit() == False)
        ;
    if (((c = bioskey(0)) & 0xff) == 0)
        c = (c >> 8) | 0x1080;
    else
        c &= 0xff;
    return c & 0x10ff;
}

/* ---------- read the keyboard shift status --------- */
int Keyboard::GetShift()
{
    regs.h.ah = 2;
    int86(KEYBRD, &regs, &regs);
    return regs.h.al;
}

/* ------ convert an Alt+ key to its letter equivalent ----- */
int Keyboard::AltConvert(int c)
{
    int i, a = 0;
    for (i = 0; i < 36; i++)
        if (c == altconvert[i])
            break;
    if (i < 26)
        a = 'a' + i;
    else if (i < 36)
        a = '0' + i - 26;
    return a;
}

Bool Keyboard::ShiftChanged()
{
    int sk = GetShift();
    Bool rtn = (Bool) (sk != shift);
    shift = sk;
    return rtn;
}

// ------ dispatch keyboard events
void Keyboard::DispatchEvent()
{
	DFWindow *Kwnd = desktop.InFocus();
	if (Kwnd == 0)
		Kwnd = (DFWindow *) desktop.ApplWnd();

    if (ShiftChanged())
           // ---- the shift status changed
        Kwnd->ShiftChanged(GetShift());
    if (KeyHit())	{
        // --- a key was pressed
        static int cvt[] = {END,DN,PGDN,BS,'5',FWD,HOME,UP,PGUP};
		int ch = GetKey();
        if (GetShift() & (LEFTSHIFT | RIGHTSHIFT))
            if (ch >= '1' && ch <= '9')
                ch = cvt[ch-'1'];
		if (ch != '\r' && (ch < ' ' || ch > 127))
			clearBIOSbuffer();
        Kwnd->Keyboard(ch);
		keydown = True;
	}
	if (keydown && (inp(0x60) & 0x80))	{
		// ----- the keyboard was released
		keydown = False;
		Kwnd->KeyReleased();
	}
}



