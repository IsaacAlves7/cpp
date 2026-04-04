// -------- border.cpp

#include "desktop.h"
#include "dfwindow.h"

void DFWindow::AdjustBorders()
{
	BorderAdj = TopBorderAdj = BottomBorderAdj = 0;
	if (Attribute() & BORDER)
		BorderAdj = TopBorderAdj = BottomBorderAdj = 1;
	if (Attribute() & TITLEBAR)
		TopBorderAdj = 1;
}

void DFWindow::Title()
{
	if (isVisible() && (Attribute() & TITLEBAR) && title != 0)	{
		int wd = ClientWidth();
		int tlen = min(title->Strlen(), wd-4);
		String sp1((wd-tlen)/2, ' ');
		String sp2(wd - (sp1.Strlen() + tlen), ' ');
		String tlin = sp1 + title->left(tlen) + sp2;
  		int fg = BLACK;

		WriteWindowString(tlin, BorderAdj, 0, fg, CYAN);
		if (Attribute() & CONTROLBOX)
			WriteWindowChar(CONTROLBOXCHAR, 2, 0, fg, CYAN);
	}
}

void DFWindow::DrawBorder(BoxLines bl)
{
	int wd = ClientWidth();
	int ht = Height() - 1;
	int rt = Width() - 1;

	String topline(wd, bl.n);
	String bottomline(wd, bl.s);

	// ------- top border
	int fg = FrameFG();
	int bg = FrameBG();
	WriteWindowChar(bl.nw, 0,  0, fg, bg);
	WriteWindowChar(bl.ne, rt, 0, fg, bg);
	if (!(Attribute() & TITLEBAR))
		WriteWindowString(topline, 1, 0, fg, bg);

	// ------ side borders
	for (int y = 1; y < ht; y++)	{
		WriteWindowChar(bl.w, 0, y, fg, bg);
		if (!(Attribute() & VSCROLLBAR))
			WriteWindowChar(bl.e, rt, y, fg, bg);
	}

	// ----- bottom border
	WriteWindowChar(bl.sw, 0,  y, fg, bg);
	WriteWindowChar(bl.se, rt, y, fg, bg);
	if (!(Attribute() & (HSCROLLBAR | STATUSBAR)))
		WriteWindowString(bottomline, 1, y, fg, bg);
}

static BoxLines inFocusBorder = {
    FOCUS_NW,
    FOCUS_LINE,
    FOCUS_NE,
    FOCUS_SIDE,
    FOCUS_SE,
    FOCUS_LINE,
    FOCUS_SW,
    FOCUS_SIDE
};

static BoxLines outFocusBorder = {
    NW,
    LINE,
    NE,
    SIDE,
    SE,
    LINE,
    SW,
    SIDE
};

void DFWindow::Border()
{
	if (isVisible())	{
		if (Attribute() & BORDER)	{
			if (dblborder && this == desktop.InFocus())
				DrawBorder(inFocusBorder);
    		else
				DrawBorder(outFocusBorder);
		}
		Title();
	}
}

void DFWindow::Shadow()
{
	if (isVisible() && (Attribute() & SHADOW))	{
		int lf = Left();
		int tp = Top();
		int ht = Height();
		int wd = Width();
		int x = wd;
		int c;
		for (int y = 1; y < ht; y++)	{
	    	c = desktop.screen().GetVideoChar(x+lf, y+tp);
			WriteWindowChar(c, x, y, ShadowFG, ShadowBG);
		}
		String ln(wd, ' ');
		for (x = 0; x < wd; x++)
	    	ln[x] = desktop.screen().GetVideoChar(x+1+lf, y+tp);
	    int lc = ln[wd-1];	// wws will truncate ln by one char
		WriteWindowString(ln, 1, y, ShadowFG, ShadowBG);
		WriteWindowChar(lc, wd, y, ShadowFG, ShadowBG);
	}
}



