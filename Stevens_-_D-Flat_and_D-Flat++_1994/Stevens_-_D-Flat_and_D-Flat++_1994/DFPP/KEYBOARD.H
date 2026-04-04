// ------------ keyboard.h

#ifndef KEYBOARD_H
#define KEYBOARD_H

#include <dos.h>
#include "dflatdef.h"

class Keyboard    {
    union REGS regs;
    int shift;		// BIOS shift state
	Bool keydown;	// True if the key is down
	Bool insert;	// insert/overstrike mode
public:
    Keyboard();
    Bool ShiftChanged();
    int AltConvert(int);
    int GetKey();
    int GetShift();
    Bool KeyHit();
    void clearBIOSbuffer();
    void DispatchEvent();
    int ShiftState() { return shift = GetShift(); }
	void SetInsertMode()   { insert = True; }
	void ClearInsertMode() { insert = False; }
	Bool InsertMode()      { return insert; }
};

const int KEYBOARDVECT = 9;
const int KEYBOARDPORT = 0x60;

inline void Keyboard::clearBIOSbuffer()
{
    *(int far *)(MK_FP(0x40,0x1a)) = 
    *(int far *)(MK_FP(0x40,0x1c));
}

// ----- keyboard BIOS (0x16) functions
const int KEYBRD = 0x16;
const int READKB = 0;
const int KBSTAT = 1;

const int ZEROFLAG = 0x40;
const int OFFSET = 0x1000;

const int RUBOUT = 8;
const int BELL   = 7;
const int ESC    = 27;
const unsigned F1          = (187+OFFSET);
const unsigned F8          = (194+OFFSET);
const unsigned SHIFT_F8    = (219+OFFSET);
const unsigned F10         = (196+OFFSET);
const unsigned HOME        = (199+OFFSET);
const unsigned UP          = (200+OFFSET);
const unsigned PGUP        = (201+OFFSET);
const unsigned BS          = (203+OFFSET);
const unsigned FWD         = (205+OFFSET);
const unsigned END         = (207+OFFSET);
const unsigned DN          = (208+OFFSET);
const unsigned PGDN        = (209+OFFSET);
const unsigned INS         = (210+OFFSET);
const unsigned DEL         = (211+OFFSET);
const unsigned CTRL_HOME   = (247+OFFSET);
const unsigned CTRL_PGUP   = (132+OFFSET);
const unsigned CTRL_BS     = (243+OFFSET);
const unsigned CTRL_FIVE   = (143+OFFSET);
const unsigned CTRL_FWD    = (244+OFFSET);
const unsigned CTRL_END    = (245+OFFSET);
const unsigned CTRL_PGDN   = (246+OFFSET);
const unsigned SHIFT_HT    = (143+OFFSET);
const unsigned ALT_A       = (158+OFFSET);
const unsigned ALT_B       = (176+OFFSET);
const unsigned ALT_C       = (174+OFFSET);
const unsigned ALT_D       = (160+OFFSET);
const unsigned ALT_E       = (146+OFFSET);
const unsigned ALT_F       = (161+OFFSET);
const unsigned ALT_G       = (162+OFFSET);
const unsigned ALT_H       = (163+OFFSET);
const unsigned ALT_I       = (151+OFFSET);
const unsigned ALT_J       = (164+OFFSET);
const unsigned ALT_K       = (165+OFFSET);
const unsigned ALT_L       = (166+OFFSET);
const unsigned ALT_M       = (178+OFFSET);
const unsigned ALT_N       = (177+OFFSET);
const unsigned ALT_O       = (152+OFFSET);
const unsigned ALT_P       = (153+OFFSET);
const unsigned ALT_Q       = (144+OFFSET);
const unsigned ALT_R       = (147+OFFSET);
const unsigned ALT_S       = (159+OFFSET);
const unsigned ALT_T       = (148+OFFSET);
const unsigned ALT_U       = (150+OFFSET);
const unsigned ALT_V       = (175+OFFSET);
const unsigned ALT_W       = (145+OFFSET);
const unsigned ALT_X       = (173+OFFSET);
const unsigned ALT_Y       = (149+OFFSET);
const unsigned ALT_Z       = (172+OFFSET);
const unsigned ALT_1      = (0xf8+OFFSET);
const unsigned ALT_2      = (0xf9+OFFSET);
const unsigned ALT_3      = (0xfa+OFFSET);
const unsigned ALT_4      = (0xfb+OFFSET);
const unsigned ALT_5      = (0xfc+OFFSET);
const unsigned ALT_6      = (0xfd+OFFSET);
const unsigned ALT_7      = (0xfe+OFFSET);
const unsigned ALT_8      = (0xff+OFFSET);
const unsigned ALT_9      = (0x80+OFFSET);
const unsigned ALT_0      = (0x81+OFFSET);
const unsigned ALT_HYPHEN = (130+OFFSET);
const unsigned CTRL_F4    = (225+OFFSET);
const unsigned ALT_F4     = (235+OFFSET);
const unsigned ALT_F6     = (237+OFFSET);

enum {CTRL_A=1,CTRL_B,CTRL_C,CTRL_D,CTRL_E,CTRL_F,CTRL_G,CTRL_H,
        CTRL_I,CTRL_J,CTRL_K,CTRL_L,CTRL_M,CTRL_N,CTRL_O,CTRL_P,
        CTRL_Q,CTRL_R,CTRL_S,CTRL_T,CTRL_U,CTRL_V,CTRL_W,CTRL_X,
        CTRL_Y,CTRL_Z };

// ------- BIOS shift key mask values
const int RIGHTSHIFT = 0x01;
const int LEFTSHIFT  = 0x02;
const int CTRLKEY    = 0x04;
const int ALTKEY     = 0x08;
const int SCROLLLOCK = 0x10;
const int NUMLOCK    = 0x20;
const int CAPSLOCK   = 0x40;
const int INSERTKEY  = 0x80;

#endif

