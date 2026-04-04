// ----------- screen.h

#ifndef SCREEN_H
#define SCREEN_H

#include <dos.h>
#include "dflatdef.h"
#include "rectangl.h"
#include "strings.h"

class Screen    {
    unsigned address;
    unsigned mode;
    unsigned page;
    unsigned height;
    unsigned width;
    union REGS regs;
    // ---- compute video offset address
    unsigned vad(int x, int y) { return y * (width*2) + x*2; }
public:
    Screen();
    unsigned Height() { return height; }
    unsigned Width()  { return width; }
    unsigned Page()   { return page; }
    Bool isEGA();
    Bool isVGA();
    Bool isMono() { return (Bool) (mode == 7); }
    Bool isText() { return (Bool) (mode < 4);  }
    void Scroll(Rect& rc, int d, int fg, int bg);
    unsigned int GetVideoChar(int x, int y);
    void PutVideoChar(int x, int y, unsigned int c);
    void WriteVideoString(const String& s,int x,int y,int fg,int bg);
    void GetBuffer(const Rect& rc, char *bf);
    void PutBuffer(const Rect& rc, char *bf);
};

const int VIDEO = 0x10;

inline int clr(int fg, int bg)
{
    return fg | (bg << 4);
}

#endif



