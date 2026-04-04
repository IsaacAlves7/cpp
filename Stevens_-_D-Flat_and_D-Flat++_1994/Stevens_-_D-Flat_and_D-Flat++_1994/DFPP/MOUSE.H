// ---------- mouse.h

#ifndef MOUSE_H
#define MOUSE_H

#include <dos.h>
#include "dfwindow.h"
#include "timer.h"

class Mouse    {
    Bool installed;      // True = mouse is installed
    char *statebuffer;   // mouse state buffer
    Timer doubletimer;   // mouse double-click timer
    Timer delaytimer;    // mouse typematic click timer
    int prevx;           // previous mouse x coordinate
    int prevy;           // previous mouse y coordinate
    int clickx;          // click x position
    int clicky;          // click y position
    int releasex;        // release x position
    int releasey;        // release y position
    union REGS regs;
    DFWindow *MouseWindow(int mx, int my);
    void CallMouse(int m1,int m2=0,int m3=0,int m4=0,unsigned es=0);
    void DispatchRelease();
    void DispatchMove();
    void DispatchLeftButton();
public:
    Mouse();
    ~Mouse();
    Bool Installed() { return installed; }
    void GetPosition(int& x, int& y); // get mouse position
    void SetPosition(int x, int y);   // set mouse position
    Bool Moved();          // True if mouse has moved
    void Show();           // show the mouse cursor
    void Hide();           // hide the mouse cursor
    Bool LeftButton();     // True if left button is pressed
    Bool ButtonReleased(); // True if button was released
    void SetTravel(int minx, int maxx, int miny, int maxy);
	void RestoreTravel();
    void DispatchEvent();
};

const int MOUSE          = 0x33;  // mouse interrupt vector
// -------- mouse commands
const int RESETMOUSE     =  0;
const int SHOWMOUSE      =  1;
const int HIDEMOUSE      =  2;
const int READMOUSE      =  3;
const int SETPOSITION    =  4;
const int BUTTONRELEASED =  6;
const int XLIMIT         =  7;
const int YLIMIT         =  8;
const int BUFFSIZE       = 21;
const int SAVESTATE      = 22;
const int RESTORESTATE   = 23;
// -------- timer delays for mouse repeat, double clicks
const int DELAYTICKS     =  1;
const int FIRSTDELAY     =  7;
const int DOUBLETICKS    =  5;

#endif


