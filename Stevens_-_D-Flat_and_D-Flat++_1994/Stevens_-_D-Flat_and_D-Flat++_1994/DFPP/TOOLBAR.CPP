// ------------ toolbar.cpp

#include "toolbar.h"

// ----- construct a Toolbar (par *must* -> a window)
ToolBar::ToolBar(DFWindow *par) :
		Control(par->ClientLeft(), par->ClientTop(),
			    2, par->ClientWidth(), par)
{
    windowtype = ToolbarWindow;
    par->SetAttribute(TOOLBAR);
	Color clr = Colors();
    clr.fg = clr.bg = ToolBarBG;
	SetColor(clr);
}

// ---- resize the menubar when the application window resizes
void ToolBar::ParentSized(int xdif, int)
{
    Size(Right()+xdif, Bottom());
}

void ToolBar::Paint()
{
	if (isEnabled())
		Control::Paint();
}

