// ------------- button.cpp

#include "button.h"
#include "desktop.h"

Button::Button(const String& lbl, int lf, int tp, DFWindow *par)
                : TextBox(lf, tp, 1, lbl.Strlen()+5, par)
{
    String lb("( ) ");
    lb += lbl;
    lb += " ";
    SetText(lb);
    setting = False;
    setchar = ' ';

	Color clr = Colors();
    clr.fg =
    clr.sfg =
    clr.ffg =
    clr.hfg = par ? par->ClientFG() : 0; 
    clr.bg =
    clr.sbg =
    clr.fbg =
    clr.hbg = par ? par->ClientBG() : 0;
	SetColor(clr);
    shortcutfg = RED;
}

void Button::Paint()
{
    if (isVisible())    {
        (*text)[1] = setting ? setchar : ' ';
        if (isEnabled())
            WriteShortcutLine(0, ClientFG(), ClientBG());
        else 
            WriteTextLine(0, HighlightFG(), HighlightBG());
    }
}

Bool Button::SetFocus()
{
    TextBox::SetFocus();
    desktop.cursor().NormalCursor();
    desktop.cursor().SetPosition(Left()+1, Top());
    desktop.cursor().Show();
    return True;
}

void Button::ResetFocus()
{
    TextBox::ResetFocus();
    desktop.cursor().Hide();
}

void Button::Keyboard(int key)
{
    if (key == ' ')
        InvertButton();
    else 
        TextBox::Keyboard(key);
}

void Button::LeftButton(int mx, int my)
{
    if (ClientRect().Inside(mx,my))
        InvertButton();
}

void Button::InvertButton()
{
    if (setting)
        ReleaseButton();
    else
        PushButton();
}

void Button::PushButton()
{
    setting = True;
    Paint();
}

void Button::ReleaseButton()
{
    setting = False;
    Paint();
}






