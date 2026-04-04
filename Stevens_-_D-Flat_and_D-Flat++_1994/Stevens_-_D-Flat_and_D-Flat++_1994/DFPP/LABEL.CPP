// ------------- label.cpp

#include <ctype.h>
#include "label.h"

void Label::InitLabel(const String& txt)
{
	SetText(String(txt));
	int off = txt.FindChar(SHORTCUTCHAR);
	if (off != -1)
		Control::LatestShortcut = tolower(txt[off+1]);
	Disable();
}

Label::Label(const String& txt, int lf, int tp, int wd, DFWindow *par)
					: TextBox(lf, tp, 1, wd ? wd : txt.Strlen(), par)
{
	InitLabel(txt);
}

Label::Label(const String& txt, int lf, int tp, DFWindow *par)
					: TextBox(lf, tp, 1, txt.Strlen(), par)
{
	InitLabel(txt);
}

// -------- paint the window 
void Label::Paint()
{
	if (Parent() != 0)
		SetColor(Parent()->Colors());
    shortcutfg = RED;
    if (isVisible())
        WriteShortcutLine(0, ClientFG(), ClientBG());
}




