// ----------- screen.cpp

#include "desktop.h"

Screen::Screen()
{
    if (isEGA() || isVGA())    {
        // --- turn blinking off
        regs.x.ax = 0x1003;
        regs.h.bl = 0;
        int86(VIDEO, &regs, &regs);
    }
    // ---- get the video mode and page
    regs.h.ah = 15;
    int86(VIDEO, &regs, &regs);
    mode = regs.h.al;
    page = regs.x.bx;
    page &= 0xff00;
    mode &= 0x7f;
    // ---- Monochrome Display Adaptor or text mode
    if (isMono())
        address = 0xb000;
    else
        // ------ Text mode
        address = 0xb800 + page;
    width = *(unsigned char far *)( MK_FP(0x40,0x4a) );
    if (isVGA() || isEGA())
        height = *(unsigned char far *)( MK_FP(0x40,0x84) )+1;
    else
        height = 25;
}

// ---- test for EGA
Bool Screen::isEGA(void)
{
    if (isVGA())
        return False;
    regs.h.ah = 0x12;
    regs.h.bl = 0x10;
    int86(VIDEO, &regs, &regs);
    return (Bool) (regs.h.bl != 0x10);
}

// ---- test for VGA
Bool Screen::isVGA(void)
{
    regs.x.ax = 0x1a00;
    int86(VIDEO, &regs, &regs);
    return (Bool) (regs.h.al == 0x1a && regs.h.bl > 6);
}

// --------- scroll the screen d: 1 = up, 0 = dn
void Screen::Scroll(Rect& rc, int d, int fg, int bg)
{
    desktop.mouse().Hide();
    regs.h.cl = rc.Left();
    regs.h.ch = rc.Top();
    regs.h.dl = rc.Right();
    regs.h.dh = rc.Bottom();
    regs.h.bh = clr(fg,bg);
    regs.h.ah = 7 - d;
    regs.h.al = 1;
    int86(VIDEO, &regs, &regs);
    desktop.mouse().Show();
}

// -------- read a character of video memory
unsigned int Screen::GetVideoChar(int x, int y)
{
    int c;
    desktop.mouse().Hide();
	int far *fp = (int far *) MK_FP(address, vad(x,y));
	c = *fp;
//    c = peek(address, vad(x,y));
    desktop.mouse().Show();
    return c & 255;
}

// -------- write a character of video memory
void Screen::PutVideoChar(int x, int y, unsigned int c)
{
    if (x >= 0 && x < width && y >= 0 && y < height)    {
        desktop.mouse().Hide();
		int far *fp = (int far *) MK_FP(address, vad(x,y));
		*fp = c;
//        poke(address, vad(x,y), c);
        desktop.mouse().Show();
    }
}

// --------- Write a string to video memory
void Screen::WriteVideoString(const String& s,int x,int y,int fg,int bg)
{
    if (x >= 0 && x < width && y >= 0 && y < height)    {
        int len = s.Strlen();
        int *ln = new int[len];
        int *cp1 = ln;
        int col = clr(fg,bg) << 8;

		for (int i = 0; i < s.Strlen(); i++)
            *cp1++ = (s[i] & 255) | col;

        if (x + len >= width)
            len = width - x;
        desktop.mouse().Hide();
        movedata(FP_SEG(ln),FP_OFF(ln),address,vad(x,y),len*2);
        desktop.mouse().Show();
        delete [] ln;
    }
}

// -- read a rectangle of video memory into a save buffer
void Screen::GetBuffer(const Rect& rc, char *bf)
{
    if (rc.Left() >= width || rc.Top() >= height)
        return;
    int ht = rc.Bottom()-rc.Top()+1;
    int bytes_row = (rc.Right()-rc.Left()+1) * 2;
    unsigned vadr = vad(rc.Left(), rc.Top());
    desktop.mouse().Hide();
    while (ht--)    {
        movedata(address, vadr, FP_SEG(bf),
                FP_OFF(bf), bytes_row);
        vadr += width*2;
        bf = (char far *)bf + bytes_row;
    }
    desktop.mouse().Show();
}

// -- write a rectangle of video memory from a save buffer
void Screen::PutBuffer(const Rect& rc, char *bf)
{
    if (rc.Left() >= width || rc.Top() >= height)
        return;
    int ht = rc.Bottom()-rc.Top()+1;
    int bytes_row = (rc.Right()-rc.Left()+1) * 2;
    unsigned vadr = vad(rc.Left(), rc.Top());
    desktop.mouse().Hide();
    while (ht--)    {
        movedata(FP_SEG(bf), FP_OFF(bf), address,
                vadr, bytes_row);
        vadr += width*2;
        bf += bytes_row;
    }
    desktop.mouse().Show();
}



