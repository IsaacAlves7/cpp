// ------------- pbutton.cpp

#include <ctype.h>
#include "pbutton.h"
#include "desktop.h"

static Color col = {
	BLACK,				// fg
	CYAN,				// bg
	WHITE,				// selected fg
	CYAN,				// selected bg
	BLACK,				// frame fg
	CYAN,				// frame bg
	DARKGRAY,			// highlighted fg
	CYAN				// highlighted bg
};

PushButton::PushButton(const String& lbl, int lf, int tp, DFWindow *par)
                : TextBox(lf, tp, 1, lbl.Strlen()+1, par)
{
    windowtype = PushButtonWindow;
	SetAttribute(SHADOW);
    pressed = False;
    cmdfunction = 0;
    String lb(" ");
    lb += lbl + " ";
    SetText(lb);
	int off = lb.FindChar(SHORTCUTCHAR);
	if (off != -1)
		shortcut = tolower(lb[off+1]);
	SetColor(col);
	shortcutfg = RED;
}


void PushButton::Paint()
{
    if (isVisible())    {
        COLORS fg;
        COLORS bg;
        if (!pressed)    {
            if (isEnabled())    {
                if (this == desktop.InFocus())
                    fg = SelectedFG();
                else
                    fg = ClientFG();
                WriteShortcutLine(0, fg, ClientBG());
            }
            else 
                WriteTextLine(0, HighlightFG(), HighlightBG());
        }
        else     {
            // ---- display a pressed button
            fg = ClientBG();
            bg = Parent()->ClientBG();
            int wd = Width();
            WriteWindowChar(' ', 0, 0, fg, bg);
            for (int x = 0; x < wd; x++)    {
                WriteWindowChar(220, x+1, 0, fg, bg);
                WriteWindowChar(223, x+1, 1, fg, bg);
            }
        }
    }
}

void PushButton::Shadow()
{
    if (isVisible() && !pressed && (Attribute() & SHADOW))    {
        COLORS bg = Parent()->ClientBG();
        int wd = Width();
        WriteWindowChar(220, wd, 0, BLACK, bg);
        for (int x = 1; x <= wd; x++)
            WriteWindowChar(223, x, 1, BLACK, bg);
    }
}

void PushButton::ResetFocus()
{
    TextBox::ResetFocus();
    Paint();
}

void PushButton::Keyboard(int key)
{
    if (key == '\r')
        PressButton();
    else 
        TextBox::Keyboard(key);
}

void PushButton::LeftButton(int mx, int my)
{
    if (ClientRect().Inside(mx,my))    {
        PressButton();
        CaptureFocus();
    }
}

void PushButton::MouseMoved(int mx, int my)
{
    if (desktop.FocusCapture() == this)
        if (ClientRect().Inside(mx,my))
            PressButton();
        else
            ReleaseButton();
}

void PushButton::ButtonReleased(int, int)
{
	ReleaseFocus();
    ButtonCommand();
}

void PushButton::KeyReleased()
{
    ButtonCommand();
}

void PushButton::ButtonCommand()
{
    if (pressed)    {
        ReleaseButton();
        if (cmdfunction != 0 && owner != 0)
            (owner->*cmdfunction)();
    }
}

void PushButton::PressButton()
{
    if (!pressed)    {
        pressed = True;
        Show();
    }
}

void PushButton::ReleaseButton()
{
    if (pressed)    {
        pressed = False;
        Show();
    }
}

void PushButton::ShortcutSelect()
{
	Control::ShortcutSelect();
	PressButton();
}





