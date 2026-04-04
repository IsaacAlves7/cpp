// ------------- mouse.cpp

#include <stdio.h>
#include "desktop.h"

// -------- mouse constructor
Mouse::Mouse()
{
    // ------- see if mouse driver is installed
	int far *ms42 = (int far *) MK_FP(0, MOUSE*4+2);
	int far *ms4  = (int far *) MK_FP(0, MOUSE*4);
    unsigned char far *ms = (unsigned char far *) MK_FP(*ms42, *ms4);
    // --- if the interrupt vector is null or points to a retf,
    //     the mouse driver is not installed
    installed = (Bool) (ms != 0 && *ms != 0xcf);

    if (installed)    {
        // --- get the mouse state buffer size
        CallMouse(BUFFSIZE);
        statebuffer = new char[regs.x.bx];
        // --- save the mouse state
        CallMouse(SAVESTATE, 0, 0,
            FP_OFF(statebuffer), FP_SEG(statebuffer));
        // --- reset the mouse
        CallMouse(RESETMOUSE);
        prevx = prevy =
        clickx = clicky =
        releasex = releasey -1;
    }
}

Mouse::~Mouse()
{
    if (installed)    {
        Hide();
        // --- restore the mouse state
        CallMouse(RESTORESTATE, 0, 0,
            FP_OFF(statebuffer), FP_SEG(statebuffer));
        delete [] statebuffer;
    }
}

void Mouse::GetPosition(int& mx, int& my)
{
    mx = my = 0;
    if (installed)    {
        CallMouse(READMOUSE);
        mx = regs.x.cx/8;
        my = regs.x.dx/8;
        if (desktop.screen().Width() == 40)
            mx /= 2;
    }
}

void Mouse::SetPosition(int x, int y)
{
    if (installed)    {
        if (desktop.screen().Width() == 40)
            x *= 2;
        CallMouse(SETPOSITION,0,x*8,y*8);
    }
}

Bool Mouse::Moved()
{
    int x, y;
    Bool rtn = False;
    if (installed)    {
        GetPosition(x, y);
        rtn = (Bool) (x != prevx || y != prevy);
        prevx = x;
        prevy = y;
    }
    return rtn;
}

void Mouse::Show()
{
    if (installed)
        CallMouse(SHOWMOUSE);
}

void Mouse::Hide()
{
    if (installed)
        CallMouse(HIDEMOUSE);
}

Bool Mouse::LeftButton()
{
    Bool rtn = False;
    if (installed)    {
        CallMouse(READMOUSE);
        rtn = (Bool) ((regs.x.bx & 1) == 1);
    }
    return rtn;
}

Bool Mouse::ButtonReleased()
{
    Bool rtn = False;
    if (installed)    {
        CallMouse(BUTTONRELEASED);
        rtn = (Bool) (regs.x.bx != 0);
    }
    return rtn;
}

void Mouse::SetTravel(int minx, int maxx, int miny, int maxy)
{
    if (installed)    {
        if (desktop.screen().Width() == 40)    {
            minx *= 2;
            maxx *= 2;
        }
        CallMouse(XLIMIT, 0, minx*8, maxx*8);
        CallMouse(YLIMIT, 0, miny*8, maxy*8);
    }
}

void Mouse::RestoreTravel()
{
	SetTravel(0, desktop.screen().Width()-1,
			  0, desktop.screen().Height()-1);
}

void Mouse::CallMouse(int m1,int m2,int m3,int m4, unsigned es)
{
    struct SREGS sregs;
    segread(&sregs);
    if (es != 0)
        sregs.es = es;
    regs.x.dx = m4;
    regs.x.cx = m3;
    regs.x.bx = m2;
    regs.x.ax = m1;
    int86x(MOUSE, &regs, &regs, &sregs);
}

// ------ get the window to send mouse events
DFWindow *Mouse::MouseWindow(int mx, int my)
{
    DFWindow *Mwnd = desktop.inWindow(mx, my);
	if (Mwnd == 0)
        Mwnd = desktop.FocusCapture();
	else if (desktop.FocusCapture() != 0)
		if (!Mwnd->isDescendedFrom(desktop.FocusCapture()))
    	    Mwnd = desktop.FocusCapture();
    return Mwnd;
}

void Mouse::DispatchRelease()
{
    if (ButtonReleased())    {
        int mx, my;
        GetPosition(mx, my);
        DFWindow *Mwnd;
        if ((Mwnd = MouseWindow(mx, my)) == 0)
            return;
        // ------- disable typematic check
        clickx = clicky = -1;
        delaytimer.DisableTimer();
        // ------- the button was released
        if (mx == releasex && my == releasey)
            // ---- same position as last left button release
            if (doubletimer.TimerRunning())
                // -- second click before double timeout
                Mwnd->DoubleClick(mx, my);
        doubletimer.SetTimer(DOUBLETICKS);
        Mwnd->ButtonReleased(mx, my);
        releasex = mx;
        releasey = my;
    }
}

void Mouse::DispatchMove()
{
    if (Moved())    {
        int mx, my;
        GetPosition(mx, my);
        DFWindow *Mwnd;
        if ((Mwnd = MouseWindow(mx, my)) != 0)    {
            Mwnd->MouseMoved(mx, my);
            clickx = clicky = -1;
        }
    }
}

void Mouse::DispatchLeftButton()
{
    if (LeftButton())    {
        int mx, my;
        GetPosition(mx, my);
        DFWindow *Mwnd;
        if ((Mwnd = MouseWindow(mx, my)) == 0)
            return;
        if (mx == clickx && my == clicky)    {
            if (delaytimer.TimedOut())    {
                // ---- button held down a while
                delaytimer.SetTimer(DELAYTICKS);
                // ---- post a typematic-like button
                Mwnd->LeftButton(mx, my);
				return;
            }
        }
		else
		    delaytimer.DisableTimer();
		if (delaytimer.TimerDisabled())    {
            // --------- new button press
   	        delaytimer.SetTimer(FIRSTDELAY);
       	    if (Mwnd->SetFocus())
           	    Mwnd->LeftButton(mx, my);
		}
        clickx = mx;
        clicky = my;
    }
}

// -------- dispatch mouse events
void Mouse::DispatchEvent()
{
    DispatchRelease();
    DispatchMove();
    DispatchLeftButton();
}



