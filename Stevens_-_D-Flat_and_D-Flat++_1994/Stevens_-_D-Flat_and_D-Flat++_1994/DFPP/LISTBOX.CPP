// ------------ listbox.cpp

#include <ctype.h>
#include "listbox.h"
#include "keyboard.h"

static Color col = {
	BLACK,			// fg
	LIGHTGRAY,		// bg
	LIGHTGRAY,		// selected fg
	BLACK,			// selected bg
	BLACK,			// frame fg
	LIGHTGRAY,		// frame bg
	LIGHTGRAY,		// inactive fg
	BLACK			// inactive bg
};

// ----------- common constructor code
void ListBox::OpenWindow()
{
    windowtype = ListboxWindow;
    selection = -1;
    addmode = False;
    anchorpoint = -1;
    selectcount = 0;
	SetColor(col);
}

void ListBox::ClearSelection()
{
    if (selection != -1)	{
		ClearTextBlock();
        WriteTextLine(selection);
		selection = -1;
	}
}

void ListBox::ClearText()
{
	TextBox::ClearText();
	OpenWindow();
}

void ListBox::SetSelection(int sel)
{
    ClearSelection();
	SelectText(0,sel,ClientWidth()-1,sel);
    if (sel >= 0 && sel < wlines)    {
        selection = sel;
		if (wtop > sel || sel >= wtop+ClientHeight())	{
			wtop = max(0, min(sel, wlines-ClientHeight()));
			Paint();
		}
		else
	        WriteTextLine(sel);
		Select();
    }
}

void ListBox::Keyboard(int key)
{
    int sel = selection; // (ClearSelection changes selection)
    switch (key)    {
        case UP:
            if (sel > 0)    {
                ClearSelection();
                if (sel == wtop)
                    ScrollDown();
                SetSelection(sel-1);
            }
            return;
        case DN:
            if (sel < wlines-1)    {
                ClearSelection();
                if (sel == wtop+ClientHeight()-1)
                    ScrollUp();
                SetSelection(sel+1);
            }
            return;
        case '\r':
            Choose();
            return;
        case FWD:
        case BS:
        case CTRL_PGUP:
        case CTRL_PGDN:
			// --- bypass these
            return;
        default:
			if (wlines && (Attribute() & SORTED) &&
					(isalpha(key) || isdigit(key)))	{
				int i = Selection();
				do	{
					if (++i == wlines)	{
						if (Selection() == -1)
							break;
						i = 0;
					}
					if (tolower(*TextLine(i)) == tolower(key))	{
						SetSelection(i);
						break;
					}
				}
				while (i != Selection());
			}
            break;
    }
    TextBox::Keyboard(key);
}

// ---------- Left mouse button was clicked
void ListBox::LeftButton(int mx, int my)
{
	if (my != prevmouseline)    {
		if (ClientRect().Inside(mx, my))    {
			int y = my - ClientTop();
            if (wlines && y < wlines-wtop)
                SetSelection(wtop+y);
        }
    }
    DFWindow::LeftButton(mx, my);
}

void ListBox::DoubleClick(int mx, int my)
{
    if (ClientRect().Inside(mx, my))    {
        my -= ClientTop();
        if (wlines && my < wlines-wtop)
            Choose();
    }
    DFWindow::DoubleClick(mx, my);
}

void ListBox::ButtonReleased(int, int)
{
    prevmouseline = -1;
}

ListBox *ListBox::This;

// --- this function is static so that qsort may take its address
int ListBox::StrCmp(const void *o1, const void *o2)
{
	const unsigned *ptraddr = This->TextPointers();
	int l1 = ((const unsigned *)o1 - ptraddr);
	int l2 = ((const unsigned *)o2 - ptraddr);
	String s1, s2;
	This->ExtractTextLine(s1, l1);
	This->ExtractTextLine(s2, l2);
	if (s1 < s2)
		return -1;
	if (s1 > s2)
		return 1;
	return 0;
}

// --------- append a line of text to the listbox buffer
void ListBox::AddText(const String& txt)
{
	TextBox::AddText(txt);
	if (wlines > 1 && (Attribute() & SORTED))	{
		// --- sort the new entry into the table
		This = this;
		unsigned *lp = TextPointers();
		unsigned *np = lp + wlines-1;
		while (lp < np)	{
			if (StrCmp(lp, np) >= 0)
				break;
			lp++;
		}
		if (lp < np)	{
			// ---- shift offsets to make room for new one
			unsigned p = *np;
			while (np > lp)	{
				*np = *(np-1);
				--np;
			}
			// --- move new offset into its ordered slot
			*lp = p;
		}
	}
}

void ListBox::BuildTextPointers()
{
	TextBox::BuildTextPointers();
	if (wlines > 1 && (Attribute() & SORTED))	{
		// --- sort the table
		This = this;
        qsort(TextPointers(), wlines, sizeof(unsigned),
			(int (*)(const void *, const void *)) &ListBox::StrCmp);
	}
}
