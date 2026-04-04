// ------------- cursor.h

#ifndef CURSOR_H
#define CURSOR_H

// ------- video BIOS (0x10) functions
const int SETCURSORTYPE = 1;
const int SETCURSOR     = 2;
const int READCURSOR    = 3;
const int HIDECURSOR    = 0x20;

const int MAXSAVES = 50;  // depth of cursor save/restore stack

class Cursor {
    // --- cursor save/restore stack
    int cursorpos[MAXSAVES];
    int cursorshape[MAXSAVES];
    int cs;             // count of cursor saves in effect
    union REGS regs;
    void Cursor::GetCursor();
public:
    Cursor();
    ~Cursor();
    void SetPosition(int x, int y);
    void GetPosition(int& x, int& y);
    void SetType(unsigned t);
    void NormalCursor() { SetType(0x0607); }
    void BoxCursor()    { SetType(0x0107); }
    void Hide();
    void Show();
    void Save();
    void Restore();
    void SwapStack();
};

inline void swap(int& a, int& b)
{
    int x = a;
    a = b;
    b = x;
}

#endif


