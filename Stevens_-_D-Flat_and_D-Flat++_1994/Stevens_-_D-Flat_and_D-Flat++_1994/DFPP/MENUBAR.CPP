// ------------- menubar.cpp

#include <ctype.h>
#include "desktop.h"
#include "menubar.h"
#include "menusel.h"

// -------- construct a menubar item
MenuBarItem::MenuBarItem(const char *Title, MenuSelection **Ms,
                         void (DFWindow::*MenuPrep)())
{
    if (Title != 0)
        title = new String(Title);
    else
        title = 0;
    ms = Ms;
    menuprep = MenuPrep;
    popdown = 0;
}

static Color col = {
	BLACK,				// fg
	LIGHTGRAY,			// bg
	BLACK,				// selected fg
	CYAN,				// selected bg
	BLACK,				// frame fg
	LIGHTGRAY,			// frame bg
	BLACK,				// highlighted fg
	LIGHTGRAY			// highlighted bg
};

// -------- construct a menubar
MenuBar::MenuBar( MenuBarItem *MenuItems, DFWindow *par) :
            TextBox( par->ClientLeft(), par->ClientTop(),
                        1, par->ClientWidth(), par)
{
    windowtype = MenubarWindow;
    menuitems = MenuItems;
    selection = -1;
    ispoppeddown = False;
    MenuBarItem *menu = menuitems;
    menucount = 0;
    oldfocus = 0;
    int off = 2;
	SetColor(col);
    shortcutfg = RED;
    while (menu->title != 0)    {
        int len = menu->title->Strlen()-1;
        menu->x1 = off;
        menu->x2 = off+len-1;
        off += len+2;
        String ttl("  ");
        ttl += *(menu->title);
        AddText(ttl);
        int n = text->Strlen()-1;
        (*text)[n] = '\0';
        menu->popdown = new PopDown(this, menu->ms);
        menu++;
        menucount++;
    }
	par->SetAttribute(MENUBAR);
}

// -------- menubar destructor
MenuBar::~MenuBar()
{
    MenuBarItem *menu = menuitems;
    while (menu->title != 0)    {
        delete menu->popdown;
        menu++;
    }
    TextBox::CloseWindow();
}

// ---- menubar gets the focus
Bool MenuBar::SetFocus()
{
    if (oldfocus == 0)
        if (desktop.InFocus() != 0)
            if (desktop.InFocus()->State() != CLOSING)
                oldfocus = desktop.InFocus();
    return TextBox::SetFocus();
}

// ---- menubar loses the focus
void MenuBar::ResetFocus()
{
    if (!ispoppeddown)    {
        SetSelection(-1);
        oldfocus = 0;
    }
    TextBox::ResetFocus();
}

// -------- paint the menubar
void MenuBar::Paint()
{
    WriteShortcutLine(0, ClientFG(), ClientBG());
    if (selection != -1)    {
        int x = menuitems[selection].x1;
        int len = menuitems[selection].x2-x+2;
        String sel = text->mid(len, x+selection);
        DisplayShortcutField(sel,x,0,SelectedFG(),SelectedBG());
    }
}

// ------- left mouse button is pressed
void MenuBar::LeftButton(int mx, int)
{
    mx -= Left();
    MenuBarItem *menu = menuitems;
    int sel = 0;
    while (menu->title != 0)    {
        if (mx >= menu->x1 && mx <= menu->x2)    {
            if (selection != sel || !ispoppeddown)    {
                if (ispoppeddown)    {
                    PopDown *pd = menuitems[selection].popdown;
                    if (pd->isOpen())
                        pd->CloseMenu(False);
                }
                Select(sel);
            }
            return;
        }
        sel++;
        menu++;
    }
    if (selection == -1)
        SetSelection(0);
}

void MenuBar::SetSelection(int sel)
{
    selection = sel;
    Paint();
}

// ----- programmed selection
void MenuBar::Select(int sel)    {
    selection = sel;
    Select();
}

// ------- user selection
void MenuBar::Select()
{
    Paint();
    ispoppeddown = True;
    MenuBarItem& mb = *(menuitems+selection);
    int lf = Left() + mb.x1;
    int tp = Top()+1;
    if (mb.menuprep != 0)	{
		DFWindow *par = Parent();
        (par->*mb.menuprep)();
	}
    mb.popdown->OpenMenu(lf, tp);
}

// ------ test for popdown accelerator key
Bool MenuBar::AcceleratorKey(int key)
{
    MenuBarItem *menu = menuitems;
    while (menu->title != 0)    {
        PopDown *pd = menu->popdown;
        if (pd->AcceleratorKey(key))
            return True;
        menu++;
    }
    return False;
}

// ------ test for menubar shortcut key
Bool MenuBar::ShortcutKey(int key)
{
    int altkey = desktop.keyboard().AltConvert(key);
    MenuBarItem *menu = menuitems;
    int sel = 0;
    while (menu->title != 0)    {
        int off = menu->title->FindChar(SHORTCUTCHAR);
        if (off != -1)    {
            String& cp = *(menu->title);
            int c = cp[off+1];
            if (tolower(c) == altkey)    {
                SetFocus();
                Select(sel);
                return True;
            }
        }
        sel++;
        menu++;
    }
    return False;
}

// -------- keystroke while menubar has the focus
void MenuBar::Keyboard(int key)
{
    if (AcceleratorKey(key))
        return;
    if (!ispoppeddown && ShortcutKey(key))
        return;
    switch (key)    {
        case F10:
            if (ispoppeddown)
                break;
            if (this != desktop.InFocus())    {
                if (selection == -1)
                    selection = 0;
                SetFocus();
                break;
            }
            // ------ fall through
        case ESC:
            ispoppeddown = False;
            SetSelection(-1);
            if (oldfocus != 0)
                oldfocus->SetFocus();
            else
                Parent()->SetFocus();
            break;
        case FWD:
            selection++;
            if (selection == menucount)
                selection = 0;
            if (ispoppeddown)
                Select();
            else
                Paint();
            break;
        case BS:
            if (selection == 0)
                selection = menucount;
            --selection;
            if (ispoppeddown)
                Select();
            else
                Paint();
            break;
        case '\r':
            if (selection != -1)
                Select();
            break;
        case ALT_F6:
            TextBox::Keyboard(key);
            break;
        default:
            break;
    }
}

// ---- resize the menubar when the application window resizes
void MenuBar::ParentSized(int xdif, int)
{
    Size(Right()+xdif, Bottom());
}


