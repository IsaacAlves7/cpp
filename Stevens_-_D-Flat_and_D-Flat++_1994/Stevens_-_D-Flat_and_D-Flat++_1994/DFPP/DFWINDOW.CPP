// ------------ dfwindow.cpp

#include "dfwindow.h"
#include "frame.h"
#include "desktop.h"

static Color col = {
	WHITE,			// fg
	BLACK,			// bg
	WHITE,			// selected fg
	BLACK,			// selected bg
	WHITE,			// frame fg
	BLACK,			// frame bg
	WHITE,			// highlighted fg
	BLACK			// highlighted bg
};

// -------- common constructor initialization code
void DFWindow::InitWindow(int lf, int tp,
                int ht, int wd, DFWindow *par)
{
    windowtype = DFlatWindow;
    if (lf == -1)
        lf = (desktop.screen().Width()-wd)/2;
    if (tp == -1)
        tp = (desktop.screen().Height()-ht)/2;
    if (ht == -1)
        ht = desktop.screen().Height();
    if (wd == -1)
        wd = desktop.screen().Width();
    attrib = restored_attrib = 0;
    title = 0;
    ctlmenu = 0;
    videosave = 0;
    visible = False;
    prevcapture = 0;
	prevfocus = 0;
    SetParent(par);
    BorderAdj = TopBorderAdj = BottomBorderAdj = 0;
    Rect rcc(lf, tp, lf+wd-1, tp+ht-1);
    restored_rc = rect = rcc;
    clearch = ' ';
    dblborder = True;
    Enqueue();
    windowstate = OPEN;
	colors = col;
}

void DFWindow::SetTitle(const String& ttl)
{
    SetAttribute(TITLEBAR);
	delete title;
    title = new String(ttl);
}

void DFWindow::InitWindow(const String& ttl, int lf, int tp,
            int ht, int wd, DFWindow *par)
{
    InitWindow(lf, tp, ht, wd, par);
	SetTitle(ttl);
}

void DFWindow::CloseWindow()
{
    windowstate = CLOSING;
    Hide();
    // ------- close window's children
    DFWindow *Wnd = First();
    while (Wnd != 0)    {
		DFWindow *Nxt = Wnd->Next();
        Wnd->CloseWindow();
        Wnd = Nxt;
    }
    // ------ delete this window's memory
    delete title;
    title = 0;
    delete videosave;
	videosave = 0;
    DeleteCtlMenu();
    if (desktop.FocusCapture() == this)
        ReleaseFocus();
    else if (this == desktop.InFocus())    {
        if (Parent() == 0)
            desktop.SetFocus(0);
		else if (Parent()->windowstate == CLOSING)	
            desktop.SetFocus(Parent());
		else 	{
       	    NextSiblingFocus();
           	if (this == desktop.InFocus())
               	if (!Parent()->SetFocus())
                   	desktop.SetFocus(0);
        }
    }
    Dequeue();
    windowstate = CLOSED;
}

// ---------- display the window
void DFWindow::Show()
{
    if (windowstate != CLOSING)	{
    	if (attrib & SAVESELF)    {
        	Rect rc = ShadowedRect();
        	if (videosave == 0)    {
            	int sz = rc.Height() * rc.Width() * 2;
            	videosave = new char[sz];
        	}
        	if (!visible)
            	desktop.screen().GetBuffer(rc, videosave);
    	}
    	visible = True;
    	Paint();
    	Border();
    	Shadow();
    	// --- show the children of this window
    	DFWindow *Wnd = First();
    	while (Wnd != 0)    {
        	Wnd->Show();
			Wnd = Wnd->Next();
    	}
	}
}

Rect DFWindow::ShadowedRect()
{
    Rect rc = rect;
    if (attrib & SHADOW)    {
        rc.SetRight(rc.Right()+1);
        rc.SetBottom(rc.Bottom()+1);
    }
    return rc;
}

void DFWindow::Hide()
{
    if (visible)    {
        visible = False;
        // ----- hide the children
        DFWindow *Wnd = First();
        while (Wnd != 0)    {
            Wnd->Hide();
            Wnd = Wnd->Next();
        }
        Rect rc = ShadowedRect();
        if (videosave != 0)    {
            desktop.screen().PutBuffer(rc, videosave);
            delete videosave;
            videosave = 0;
        }
    }
}

Bool DFWindow::isDescendedFrom(DFWindow *ancestor)
{
	DFWindow *wnd = this;
	while (wnd != 0)	{
		if (wnd == ancestor)
			return True;
		wnd = wnd->Parent();
	}
	return False;
}

void DFWindow::Keyboard(int key)
{
    switch (key)    {
		case F1:
			HelpFunction();
			return;
		case ALT_F6:
			NextSiblingFocus();
			return;
		case ' ':
			if ((desktop.keyboard().GetShift() & ALTKEY) == 0)
				break;
			// --- fall through
        case ALT_HYPHEN:
			if (attrib & CONTROLBOX)	{
	            OpenCtlMenu();
    	        return;
			}
			break;
		case CTRL_F4:
			if (Parent() && (attrib & CONTROLBOX))	{
				CloseWindow();
    	        return;
			}
			break;
        default:
            break;
    }
    // --- send all unprocessed keystrokes 
    //     to the parent window
    if (Parent() != 0)
        Parent()->Keyboard(key);
}

void DFWindow::HelpFunction()
{
	if (Parent() != 0)
		Parent()->HelpFunction();
}

void DFWindow::ShiftChanged(int sk)
{
    if (Parent() != 0)
        Parent()->ShiftChanged(sk);
}

void DFWindow::DoubleClick(int mx, int my)
{
    if (HitControlBox(mx, my))
        CloseWindow();
}

void DFWindow::LeftButton(int mx, int my)
{
    if (my == Top())    {
        // ----- hit the top border
        if (HitControlBox(mx, my) && (attrib & CONTROLBOX))
            // ------- hit the control box
            OpenCtlMenu();
        else if (attrib & MOVEABLE)
            // ---- move the window
            desktop.frame().OpenFrame(this, True, mx);
    }
    else if ((attrib & SIZEABLE) && windowstate == OPEN)
        if (mx == Right() && my == Bottom())
            // --- hit the lower right corner, size the window
            desktop.frame().OpenFrame(this, False, 0);
    prevmouseline = my;
    prevmousecol = mx;
}

void DFWindow::ButtonReleased(int, int)
{
    prevmouseline = -1;
    prevmousecol = -1;
}

Rect DFWindow::ClientRect()
{
    Rect rc(ClientLeft(), ClientTop(),
        ClientRight(), ClientBottom());
    return rc;
}

void DFWindow::ChangePosition(int x, int y)
{
    int ht = Height();
    int wd = Width();
    rect.SetLeft(x);
    rect.SetTop(y);
    rect.SetRight(Left()+wd-1);
    rect.SetBottom(Top()+ht-1);
    restored_rc = rect;
}

// ------------ move a window
void DFWindow::Move(int x, int y)
{
    int xdif = x - Left();
    int ydif = y - Top();
    if (xdif == 0 && ydif == 0)
        return;
    Bool wasVisible = visible;
    if (wasVisible)
        Hide();
	ChangePosition(x, y);
    DFWindow *Wnd = First();
    while (Wnd != 0)    {
        Wnd->Move(Wnd->Left()+xdif, Wnd->Top()+ydif);
        Wnd = Wnd->Next();
    }
    if (wasVisible)
        Show();
}

// ------------ size a window
void DFWindow::Size(int x, int y)
{
    int xdif = x - Right();
    int ydif = y - Bottom();
    if (xdif == 0 && ydif == 0)
        return;
    Bool wasVisible = visible;
    if (wasVisible)
        Hide();
    rect.SetRight(x);
    rect.SetBottom(y);
    restored_rc = rect;

    DFWindow *Wnd = First();
    while (Wnd != 0)    {
        Wnd->ParentSized(xdif, ydif);
        Wnd = Wnd->Next();
    }
    if (wasVisible)
        Show();
}

