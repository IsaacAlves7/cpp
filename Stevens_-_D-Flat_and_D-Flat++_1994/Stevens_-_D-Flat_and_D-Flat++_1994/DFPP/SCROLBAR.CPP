// ------------- scrolbar.cpp

#include "desktop.h"
#include "textbox.h"
#include "scrolbar.h"

ScrollBar::ScrollBar(BarPlane Plane, TextBox *par)
						: DFWindow(1, 1, par)
{
	windowtype = ScrollbarWindow;
	if (par == 0)
		return;
	scrollbox = 1;
	sliding = False;
	SetClearChar(SCROLLBARCHAR);
	plane = Plane;
	if (plane == HORIZONTAL)	{
		Move(par->Left()+1, par->Bottom());
		Size(par->Right()-1, par->Bottom());
	}
	else	{
		Move(par->Right(), par->Top()+1);
		Size(par->Right(), par->Bottom()-1);
	}
	Color clr = Colors();
	clr.fg = par ? par->FrameFG() : 0;
	clr.bg = par ? par->FrameBG() : 0;
	SetColor(clr);
}

Bool ScrollBar::SetFocus()
{
	if (Parent() != 0)
		Parent()->SetFocus();
	return True;
}

void ScrollBar::Paint()
{
	if (isVisible())	{
		int fg = ClientFG();
		int bg = ClientBG();
		DFWindow::Paint();
		if (plane == HORIZONTAL)	{
			WriteWindowChar(LEFTSCROLLBOX, 0, 0, fg, bg);
			WriteWindowChar(RIGHTSCROLLBOX, Width()-1, 0, fg, bg);
			WriteWindowChar(SCROLLBOXCHAR, scrollbox, 0, fg, bg);
		}
		else	{
			WriteWindowChar(UPSCROLLBOX, 0, 0, fg, bg);
			WriteWindowChar(DOWNSCROLLBOX, 0, Height()-1, fg, bg);
			WriteWindowChar(SCROLLBOXCHAR, 0, scrollbox, fg, bg);
		}
	}
}

void ScrollBar::LeftButton(int mx, int my)
{
	if (Parent() != 0 && !sliding)	{
		desktop.mouse().SetTravel(Left(),
									Right(),
									Top(),
									Bottom());
		TextBox& Par = *(TextBox *)Parent();
		if (plane == VERTICAL)	{
			// -- test for hitting the vertical scroll buttons
			if (my == WindowRect().Top())
				Par.ScrollDown();
			else if (my == WindowRect().Bottom())
				Par.ScrollUp();
			// ------- test for hitting the vertical scroll box
			else if (my-WindowRect().Top() == scrollbox)	{
				sliding = True;
				desktop.mouse().SetTravel(WindowRect().Left(), WindowRect().Right(),
								WindowRect().Top()+1, WindowRect().Bottom()-1);
			}
			else	{
				// ----- hit in the scroll bar
				if (my-WindowRect().Top() < scrollbox)
					Par.PageUp();
				else
					Par.PageDown();
			}
		}
		else	{
			// -- test for hitting the horizontal scroll buttons
			if (mx == WindowRect().Left())
				Par.ScrollRight();
			else if (mx == WindowRect().Right())
				Par.ScrollLeft();
			// ------- test for hitting the horizontal scroll box
			else if (mx-WindowRect().Left() == scrollbox)	{
				sliding = True;
				desktop.mouse().SetTravel(WindowRect().Left()+1, WindowRect().Right()-1,
								WindowRect().Top(), WindowRect().Bottom());
			}
			else	{
				// ----- hit in the scroll bar
				if (mx-WindowRect().Left() < scrollbox)
					Par.PageLeft();
				else
					Par.PageRight();
			}
		}
	}
}

void ScrollBar::MouseMoved(int mx, int my)
{
	if (sliding)	{
		if (plane == VERTICAL)
			MoveScrollBox(my - WindowRect().Top());
		else
			MoveScrollBox(mx - WindowRect().Left());
	}
}

void ScrollBar::ButtonReleased(int, int)
{
	desktop.mouse().RestoreTravel();
	if (sliding && Parent() != 0)	{
		TextBox& Par = *(TextBox *)Parent();
		int pct = (scrollbox-1)*100;
		if (plane == VERTICAL)	{
			pct /= (WindowRect().Height()-2);
			Par.VerticalPagePosition(pct);
		}
		else	{
			pct /= (WindowRect().Width()-2);
			Par.HorizontalPagePosition(pct);
		}
	}
	sliding = False;
}

void ScrollBar::ParentSized(int xdif, int ydif)
{
	if (plane == HORIZONTAL)	{
		Size(Right()+xdif, Top());
		Move(Left(), Bottom()+ydif);
	}
	else	{
		Size(Left(), Bottom()+ydif);
		Move(Right()+xdif, Top());
	}
}

void ScrollBar::TextPosition(int TxPct)
{
	int len = plane == HORIZONTAL ? Width()-2 : Height()-2;
	int sb = 1 + ((len * TxPct) / 100);
	MoveScrollBox(sb);
}

void ScrollBar::MoveScrollBox(int sb)
{
	if (sb != scrollbox)	{
		int fg = ClientFG();
		int bg = ClientBG();
		int x1 = 0, y1 = 0, x2 = 0, y2 = 0;
		if (plane == HORIZONTAL)
			x1 = sb, x2 = scrollbox;
		else
			y1 = sb, y2 = scrollbox;
		WriteWindowChar(SCROLLBOXCHAR, x1, y1, fg, bg);
		WriteWindowChar(SCROLLBARCHAR, x2, y2, fg, bg);
		scrollbox = sb;
	}
}

