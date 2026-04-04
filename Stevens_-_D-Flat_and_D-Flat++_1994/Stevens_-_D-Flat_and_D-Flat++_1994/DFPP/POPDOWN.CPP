// ------------- popdown.cpp

#include <ctype.h>
#include "applicat.h"
#include "popdown.h"
#include "menusel.h"

static Color col = {
	BLACK,			// fg
	CYAN,			// bg
	BLACK,			// selected fg
	LIGHTGRAY,		// selected bg
	BLACK,			// frame fg
	CYAN,			// frame bg
	DARKGRAY,		// inactive fg
	CYAN			// inactive bg
};

PopDown::PopDown(DFWindow *par, MenuSelection **Selections)
                        : ListBox(5, 5, par)
{
    windowtype = PopdownWindow;
	selections = Selections;
    SetAttribute(BORDER | SHADOW | SAVESELF);
    selection = 0;
	SetSingleBorder();
    isopen = False;
    iscascaded = False;
    if (selections != 0)    {
        MenuDimensions();
        for (int i = 0; i < menuheight; i++)    {
            MenuSelection& ms = **(selections+i);
            BuildMenuLine(i);
            if (ms.type == CASCADER)    {
                ms.cascade = new PopDown(this, ms.cascaders);
                ms.cascade->iscascaded = True;
            }
        }
        Size(Left() + menuwidth, Top() + menuheight + 1);
    }
	SetColor(col);
    shortcutfg = RED;
}

// ---- shut down a popdown menu
void PopDown::CloseWindow()
{
    if (selections != 0)    {
        // --- delete all cascader popdowns
        for (int i = 0; selections[i]; i++)    {
            MenuSelection& ms = *selections[i];
            if (ms.type == CASCADER && ms.cascade != 0)
                delete ms.cascade;
        }
    }
    ListBox::CloseWindow();
}

// ------- pop down the menu
void PopDown::OpenMenu(int left, int top)
{
    Rect rc(0, 0, desktop.screen().Width()-1,
        desktop.screen().Height()-1);
    DFWindow *Wnd = Parent();
    while (Wnd != 0 && Wnd->WindowType() == PopdownWindow)
        Wnd = Wnd->Parent();
    if (Wnd != 0 && (Wnd = Wnd->Parent()) != 0)    {
        Rect rc = Wnd->ClientRect();
		if (Parent()->WindowType() == MenubarWindow)
			if (Wnd->WindowType() == ApplicationWindow)
				if (Wnd->Attribute() & TOOLBAR)
					rc = Rect(rc.Left(), rc.Top()-3, rc.Right(), rc.Bottom());
        left = min(max(left, rc.Left()),
                    rc.Right() - ClientWidth());
        top = min(max(top, rc.Top()),
                    rc.Bottom() - ClientHeight());
    }
    left = min(max(left, rc.Left()),
                    rc.Right()-ClientWidth()-1);
    top = min(max(top, rc.Top()),
                    rc.Bottom()-ClientHeight()-1);
    isopen = True;
    Move(left, top);
    CaptureFocus();
}

// ---------- deactivate the popdown menu
void PopDown::CloseMenu(Bool SendESC)
{
    if (isopen)    {
        // ------- close any open cascaded menus
        PopDown *Wnd = (PopDown *)First();
        while (Wnd != 0)    {
            Wnd->CloseMenu();
            Wnd = (PopDown *) (Wnd->Next());
        }
        Hide();
        isopen = False;
        ReleaseFocus();
        if (Parent() && !iscascaded && SendESC)
            Parent()->Keyboard(ESC);
    }
}

void PopDown::Show()
{
    if (isopen)
        ListBox::Show();
}

// -------- build a menu line
void PopDown::BuildMenuLine(int sel)
{
    int wd = menuwidth;
    String *ln;
    if (selections[sel]->type == SEPARATOR)
        ln = new String(--wd, LINE);
    else    {
        ln = new String(" ");
        *ln += selections[sel]->label;
        int r = wd-ln->Strlen();
        *ln += String(r, ' ');
        if (selections[sel]->type == CASCADER)
            (*ln)[wd-1] = CASCADEPOINTER;
    }
    AddText(*ln);
	delete ln;
}

// -------- compute menu width
void PopDown::MenuDimensions()
{
    int txlen = 0;
    for (int i = 0; selections[i] != 0; i++)  {
        if (selections[i]->type != SEPARATOR)    {
            int lblen = selections[i]->label.Strlen()-1;
            txlen = max(txlen, lblen);
        }
    }
    menuwidth = txlen+4;
    menuheight = i;
}

// ------ display a menu line
void PopDown::DisplayMenuLine(int lno)
{
    if (isopen)    {
        int fg, bg;
        int isActive = selections[lno]->isEnabled();
        int sfg = shortcutfg;
        if (lno == selection)    {
            fg = SelectedFG();
            bg = SelectedBG();
        }
        else if (isActive)    {
            fg = ClientFG();
            bg = ClientBG();
        }
        else     {
            fg = HighlightFG();
            bg = HighlightBG();
        }
        if (!isActive)
            shortcutfg = fg;
        WriteShortcutLine(lno, fg, bg);
        shortcutfg = sfg;
    }
}

// ------ set no selection current
void PopDown::ClearSelection()
{
    if (selection != -1)    {
        int sel = selection;
        selection = -1;
        DisplayMenuLine(sel);
    }
}

// ------ set a current menu selection
void PopDown::SetSelection(int sel)
{
    ClearSelection();
    if (sel >= 0 && sel < wlines)    {
        selection = sel;
        DisplayMenuLine(sel);
    }
}

// ---------- paint the menu
void PopDown::Paint()
{
    if (text == 0)
        ListBox::Paint();
    else    {
        for (int i = 0; i < wlines; i++)    {
            if (selections[i]->type == TOGGLE)    {
                int off = TextLine(i) - (const char *)*text;
                if (selections[i]->toggle == On)
                    (*text)[off] = CheckMark();
                else
                    (*text)[off] = ' ';
            }
            DisplayMenuLine(i);
        }
    }
}

// --------- paint the menu's border
void PopDown::Border()
{
    if (isopen && isVisible())    {
        int fg = FrameFG();
        int bg = FrameBG();
        int rt = Width()-1;
        ListBox::Border();
        for (int i = 0; i < wlines; i++)    {
            if (selections[i]->type == SEPARATOR)    {
                WriteWindowChar(LEDGE, 0, i+1, fg, bg);
                WriteWindowChar(REDGE, rt, i+1, fg, bg);
            }
        }
    }
}

// ------- test for a menu selection accelerator key
Bool PopDown::AcceleratorKey(int key)
{
    for (int i = 0; i < wlines; i++)    {
        MenuSelection& ms = **(selections+i);
        if (key == ms.accelerator)    {
            SetSelection(i);
            Choose();
            return True;
        }
    }
    return False;
}

// ------- test for a menu selection shortcut key
Bool PopDown::ShortcutKey(int key)
{
    key = tolower(key);
    for (int i = 0; i < wlines; i++)    {
        MenuSelection& ms = **(selections+i);
        int off = ms.label.FindChar(SHORTCUTCHAR);
        if (off != -1)    {
            String& cp = ms.label;
            int c = cp[off+1];
            if (key == tolower(c))    {
                SetSelection(i);
                Choose();
                return True;
            }
        }
    }
    return False;
}

// ----- keystroke while menu is popped down
void PopDown::Keyboard(int key)
{
    if (AcceleratorKey(key))
        return;
    if (ShortcutKey(key))
        return;
    switch (key)    {
        case UP:
            if (selection == 0)    {
                SetSelection(wlines-1);
                return;
            }
            if (selections[selection-1]->type == SEPARATOR)  {
                SetSelection(selection-2);
                return;
            }
            break;
        case DN:
            if (selection == wlines-1)    {
                SetSelection(0);
                return;
            }
            if (selections[selection+1]->type == SEPARATOR)  {
                SetSelection(selection+2);
                return;
            }
            break;
        case ESC:
            CloseMenu(ParentisMenu());
            return;
        case FWD:
        case BS:
            CloseMenu();
            if (Parent() != 0)    {
                Parent()->Keyboard(key);
                return;
            }
            break;
        default:
            break;
    }
    ListBox::Keyboard(key);
}

// ----- shift key status changed
void PopDown::ShiftChanged(int sk)
{
    if (sk & ALTKEY)
        CloseMenu(ParentisMenu());
}

// ---------- Left mouse button was clicked
void PopDown::LeftButton(int mx, int my)
{
    if (ClientRect().Inside(mx, my))    {
        if (my != prevmouseline)    {
            int y = my - ClientTop();
            if (selections[y]->type != SEPARATOR)
                SetSelection(y);
        }
    }
    else if (!WindowRect().Inside(mx, my))    {
        if (Parent() && my == Parent()->Bottom())
            Parent()->LeftButton(mx, my);
    }
    prevmouseline = my;
    prevmousecol = mx;
}

// ---------- Left mouse button was double-clicked
void PopDown::DoubleClick(int mx, int my)
{
    if (!WindowRect().Inside(mx, my))    {
        CloseMenu();
        if (Parent())
            Parent()->DoubleClick(mx, my);
    }
}

// ---------- Left mouse button was released
void PopDown::ButtonReleased(int mx, int my)
{
    if (ClientRect().Inside(mx, my))    {
        if (prevmouseline == my && prevmousecol == mx)
            if (selections[my-ClientTop()]->type != SEPARATOR)
                Choose();
    }
    else if (!WindowRect().Inside(mx, my))    {
        DFWindow *Wnd = desktop.inWindow(mx, my);
        if (!(Wnd == Parent() && my == Top()-1 &&
                mx >= Left() && mx <= Right()))    {
            CloseMenu(ParentisMenu());
            if (Wnd != 0 && Wnd != desktop.InFocus())
                Wnd->SetFocus();
        }
    }
}

// --------- user chose a menu selection
void PopDown::Choose()
{
    MenuSelection& ms = *selections[selection];
    if (ms.isEnabled())    {
        if (ms.type == CASCADER && ms.cascade != 0)
            // -------- cascaded menu
            ms.cascade->OpenMenu(Right(), Top()+selection);
        else    {
            if (ms.type == TOGGLE)    {
                // ---- toggle selection
                ms.InvertToggle();
                const char *cp = TextLine(selection);
				int off = cp - (const char *)*text;
                if (*cp == CheckMark())
                    (*text)[off] = ' ';
                else
                    (*text)[off] = CheckMark();
                DisplayMenuLine(selection);
            }
            if (ms.cmdfunction != 0)    {
                // ---- there is a function associated
                DFWindow *wnd = (PopDown *) this;
                // --- close all menus
                while (wnd &&
                        wnd->WindowType() == PopdownWindow) {
                    ((PopDown *)wnd)->CloseMenu();
                    wnd = (wnd->Parent());
                }
                if (wnd && wnd->WindowType() == MenubarWindow){
					if (wnd == desktop.InFocus())
	                    wnd->Keyboard(ESC);
                    wnd = (wnd->Parent());
                }
                if (wnd)
                    // ---- execute the function
                    (((Application *)wnd)->*ms.cmdfunction)();
            }
        }
    }
    else
        desktop.speaker().Beep();    // disabled selection
}

inline Bool isMenu(DFWindow *wnd)
{
    if (wnd != 0)    {
        WndType wt = wnd->WindowType();
        return (Bool) (wt==MenubarWindow || wt==PopdownWindow);
    }
    return False;
}

// ----- test for the parent as menu or menubar
Bool PopDown::ParentisMenu(DFWindow& wnd)
{
    return isMenu(wnd.Parent());
}

