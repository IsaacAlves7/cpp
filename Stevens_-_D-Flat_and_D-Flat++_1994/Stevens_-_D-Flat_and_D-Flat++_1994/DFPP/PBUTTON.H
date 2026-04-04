// ----------- pbutton.h

#ifndef PBUTTON_H
#define PBUTTON_H

#include "textbox.h"

class PushButton : public TextBox    {
    DFWindow *owner;  // window that gets the command
    void (DFWindow::*cmdfunction)();    // selection function
protected:
    Bool pressed;
public:
    PushButton(const String& lbl, int lf, int tp, DFWindow *par=0);
    // -------- API messages
    virtual void ResetFocus();
    virtual void Paint();
    virtual void Shadow();
    virtual void Keyboard(int key);
    virtual void LeftButton(int mx, int my);
    virtual void ButtonReleased(int mx, int my);
    virtual void MouseMoved(int mx, int my);
    virtual void KeyReleased();
    virtual void PressButton();
    virtual void ReleaseButton();
    virtual void ButtonCommand();
	virtual void ShortcutSelect();
    void SetButtonFunction(DFWindow *wnd,
                    void (DFWindow::*cmdfunc)())
        { owner = wnd; cmdfunction = cmdfunc; }
};

class OKButton : public PushButton	{
public:
	OKButton(int x, int y, DFWindow *par=0) :
					PushButton("  ~Ok  ", x, y, par)
		{ SetButtonFunction(par,
				&DFWindow::OKFunction); }
};

class CancelButton : public PushButton	{
public:
	CancelButton(int x, int y, DFWindow *par=0) :
					PushButton("~Cancel", x, y, par)
		{ SetButtonFunction(par, &DFWindow::CancelFunction); }
};

class YesButton : public PushButton	{
public:
	YesButton(int x, int y, DFWindow *par=0) :
					PushButton(" ~Yes ", x, y, par)
		{ SetButtonFunction(par, &DFWindow::OKFunction); }
};

class NoButton : public PushButton	{
public:
	NoButton(int x, int y, DFWindow *par=0) :
					PushButton(" ~No  ", x, y, par)
		{ SetButtonFunction(par, &DFWindow::CancelFunction); }
};

class HelpButton : public PushButton	{
public:
	HelpButton(int x, int y, DFWindow *par=0) :
					PushButton(" ~Help ", x, y, par)
		{ SetButtonFunction(par, &DFWindow::HelpFunction); }
};

#endif


