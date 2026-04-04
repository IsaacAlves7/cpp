// ------------ cursor.cpp

#include <dos.h>
#include "cursor.h"
#include "desktop.h"

Cursor::Cursor()
{
    cs = 0;
    Save();
}

Cursor::~Cursor()
{
    Restore();
}

// ------ get cursor shape and position
void Cursor::GetCursor()
{
    regs.h.ah = READCURSOR;
    regs.x.bx = desktop.screen().Page();
    int86(VIDEO, &regs, &regs);
}

// -------- get the current cursor position
void Cursor::GetPosition(int& x, int& y)
{
    GetCursor();
    x = regs.h.dl;
    y = regs.h.dh;
}

// ------ position the cursor
void Cursor::SetPosition(int x, int y)
{
    regs.x.dx = ((y << 8) & 0xff00) + x;
    regs.h.ah = SETCURSOR;
    regs.x.bx = desktop.screen().Page();
    int86(VIDEO, &regs, &regs);
}

// ------ save the current cursor configuration
void Cursor::Save()
{
    if (cs < MAXSAVES)    {
        GetCursor();
        cursorshape[cs] = regs.x.cx;
        cursorpos[cs] = regs.x.dx;
        cs++;
    }
}

// ---- restore the saved cursor configuration
void Cursor::Restore()
{
    if (cs)    {
        --cs;
        regs.x.dx = cursorpos[cs];
        regs.h.ah = SETCURSOR;
        regs.x.bx = desktop.screen().Page();
        int86(VIDEO, &regs, &regs);
        SetType(cursorshape[cs]);
    }
}

/* ---- set the cursor type ---- */
void Cursor::SetType(unsigned t)
{
    regs.h.ah = SETCURSORTYPE;
    regs.x.bx = desktop.screen().Page();
    regs.x.cx = t;
    int86(VIDEO, &regs, &regs);
}

/* ----- swap the cursor stack ------- */
void Cursor::SwapStack()
{
    if (cs > 1)    {
        swap(cursorpos[cs-2], cursorpos[cs-1]);
        swap(cursorshape[cs-2], cursorshape[cs-1]);
    }
}

/* ------ hide the cursor ------ */
void Cursor::Hide()
{
    GetCursor();
    regs.h.ch |= HIDECURSOR;
    regs.h.ah = SETCURSORTYPE;
    int86(VIDEO, &regs, &regs);
}

/* ------ show the cursor ------ */
void Cursor::Show()
{
    GetCursor();
    regs.h.ch &= ~HIDECURSOR;
    regs.h.ah = SETCURSORTYPE;
    int86(VIDEO, &regs, &regs);
}



