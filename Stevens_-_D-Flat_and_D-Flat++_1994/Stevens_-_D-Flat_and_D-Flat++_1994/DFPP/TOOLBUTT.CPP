// ---------- toolbutt.cpp

#include "toolbar.h"
#include "desktop.h"

static Bool settingfocus;

// ------- various color patterns
static Color EnabledColor = {
    BLACK,            // fg
    ToolButtonBG      // bg
};

static Color DisabledColor = {
    LIGHTGRAY,        // fg
    ToolButtonBG      // bg
};

// --- common constructor code
void ToolButton::InitWindow(const String& lbl)
{
    oldFocus = 0;
    windowtype = ToolButtonWindow;
    Size(Left()+ButtWidth, Top()+1);
    ClearAttribute(SHADOW);
	SetAttribute(SAVESELF);
	String lb(" ");
	label = lb + lbl;
    SetColors();
}

// --------- construct a Toolbar button specifying position
ToolButton::ToolButton(const String& lbl, int lf, int tp, DFWindow *par)
                : PushButton(lbl, lf, tp, par)
{
    InitWindow(lbl);
	lf = Left();
	tp = Top();
}

// ------ construct a Toolbar button, self-positioning
ToolButton::ToolButton(const String& lbl, DFWindow *par)
                : PushButton(lbl, 0, 0, par)
{
    InitWindow(lbl);
    if (par != 0 && par->WindowType() == ToolbarWindow)    {
        int btcount = 0;
        DFWindow *Wnd = par->First();
        while (Wnd != 0 && Wnd != this)    {
            if (Wnd->WindowType() == ToolButtonWindow)
                btcount++;
            Wnd = Wnd->Next();
        }
        Move(par->Left()+1+btcount*(ButtWidth + 2), par->Top());
    }
	lf = Left();
	tp = Top();
}

void ToolButton::ResetFocus()
{
	// --- avoid the PushButton class's extra Paint
    TextBox::ResetFocus();
}

// -------- set the button's colors
void ToolButton::SetColors()
{
    SetColor(isEnabled() ? EnabledColor : DisabledColor);
}

// ------- paint the button
void ToolButton::Paint()
{
	ToolBar *par = (ToolBar *) Parent();
	if (!(par->isEnabled()))
		return;
    SetColors();
	if (pressed)	{
		String lb("\n ");
		lb += label;
		SetText(lb);
		TextBox::Paint();
		String bt(ButtWidth, PRESSED_NORTH);
		WriteClientString(bt, 1, 0, DARKGRAY, ClientBG());
		WriteChar(PRESSED_WEST, Left(), Top(),
					DARKGRAY, ToolBarBG);
		WriteChar(PRESSED_WEST, Left(), Top()+1,
					DARKGRAY, ToolBarBG);
	}
	else if (!settingfocus)	{
		SetText(label);
		TextBox::Paint();
		String bt(ButtWidth, UNPRESSED_SOUTH);
		WriteClientString(bt, 0, 1, DARKGRAY, ClientBG());
		WriteChar(UNPRESSED_EAST, Left()+ButtWidth, Top(),
					DARKGRAY, ToolBarBG);
		WriteChar(UNPRESSED_EAST, Left()+ButtWidth, Top()+1,
					DARKGRAY, ToolBarBG);
	}
}

// -------- the button was pressed
void ToolButton::ButtonCommand()
{
    PushButton::ButtonCommand();
    if (oldFocus != 0)
        oldFocus->SetFocus();
	oldFocus = 0;
}

// ---- remember who had the focus before the button got it
Bool ToolButton::SetFocus()
{
    if (oldFocus == 0)
        oldFocus = desktop.InFocus();
	settingfocus = True;
    Bool rtn = PushButton::SetFocus();
	settingfocus = False;
	return rtn;
}

