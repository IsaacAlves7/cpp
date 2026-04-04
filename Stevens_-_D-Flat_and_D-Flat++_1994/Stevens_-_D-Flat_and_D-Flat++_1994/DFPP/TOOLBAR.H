// -------- toolbar.h

#ifndef TOOLBAR_H
#define TOOLBAR_H

#include "pbutton.h"

const COLORS ToolBarBG = BLUE;
const COLORS ToolButtonBG = CYAN;

// ------ Toolbar class
class ToolBar : public Control {
    void ParentSized(int xdif, int ydif);
    virtual void Paint();
	virtual Bool SetFocus() { return False; }
public:
    ToolBar(DFWindow *par);
};

// ---- pressed and unpressed frame characters
const int PRESSED_WEST   = 222;
const int PRESSED_NORTH   = 223;
const int UNPRESSED_EAST = 221;
const int UNPRESSED_SOUTH = 220;

const int ButtWidth = 6;

// ------- Toolbar button class
class ToolButton : public PushButton    {
    void SetColors();
    void InitWindow(const String& lbl);
	String label;
    DFWindow *oldFocus;
    virtual void ButtonCommand();
    virtual Bool SetFocus();
	int lf, tp;
public:
    ToolButton(const String& lbl, int lf, int tp, DFWindow *par=0);
    ToolButton(const String& lbl, DFWindow *par=0);
    virtual void Paint();
    virtual void ResetFocus();
};

#endif
