// ---------- paint.cpp

#include "desktop.h"
#include "dfwindow.h"

void DFWindow::ClearWindow()
{
	String ln(ClientWidth(), clearch);

	int h = ClientHeight();
	for (int y = 0; y < h; y++)
		WriteClientString(ln, 0, y, colors.fg, colors.bg);
}

void DFWindow::Paint()
{
	if (visible)
		ClearWindow();
}

void DFWindow::WriteString(const String& ln, int x, int y, int fg, int bg)
{
    int xx = x - Left();
	String line(ln);
    if (line.Strlen()+xx > Width())
        line.ChangeLength(Width()-xx);
	if (this->isDescendedFrom(desktop.InFocus()))
		desktop.screen().WriteVideoString(line, x, y, fg, bg);
	else 	{
		const char *cp = line;
		while (*cp)
			WriteChar(*cp++, x++, y, fg, bg);
	}
}

void DFWindow::WriteClientString(const String& ln,
                    int x,int y,int fg,int bg)
{
	String line(ln);
	if (line.Strlen()+x > ClientWidth())
        line.ChangeLength(ClientWidth()-x);
    WriteString(line,x+ClientLeft(),y+ClientTop(),fg,bg);
}

void DFWindow::WriteChar(int ch, int x, int y, int fg, int bg)
{
	if (ShadowedRect().Inside(x, y))	{
	    // --- coordinates are inside target window
		// --- see if overlapped by sibling
		DFWindow *sibl = Next();
		while (sibl != 0)	{
			if (sibl->isVisible() && sibl->ShadowedRect().Inside(x, y))
				return;
			sibl = sibl->Next();
		}
	   	desktop.screen().PutVideoChar(x, y, (ch & 255) | (clr(fg,bg) << 8));
	}
}


