// ------------- editbox.cpp

#include <ctype.h>
#include "desktop.h"
#include "editbox.h"

// ----------- common constructor code
void EditBox::OpenWindow()
{
	windowtype = EditboxWindow;
	column = 0;
	changed = False;
	text = new String(1,0);
	BuildTextPointers();
    insertmode = desktop.keyboard().InsertMode();
}

Bool EditBox::SetFocus()
{
	Bool rtn = TextBox::SetFocus();
	if (rtn)
		ResetCursor();
	return rtn;
}

void EditBox::ResetFocus()
{
	desktop.cursor().Hide();
	TextBox::ResetFocus();
}

// -------- process keystrokes
void EditBox::Keyboard(int key)
{
	int shift = desktop.keyboard().GetShift();
	if ((shift & ALTKEY) == 0)	{
    	switch (key)    {
        	case HOME:
            	Home();
				TestMarking();
				return;
        	case END:
            	End();
				TestMarking();
				return;
        	case CTRL_FWD:
            	NextWord();
				TestMarking();
				return;
        	case CTRL_BS:
            	PrevWord();
				TestMarking();
				return;
        	case FWD:
            	Forward();
				TestMarking();
				return;
        	case BS:
            	Backward();
				TestMarking();
				return;
        	case RUBOUT:
				ClearTextBlock();
				if (CurrentCharPosition() == 0)
					break;
            	Backward();
				// --- fall through
        	case DEL:
				if (key != RUBOUT && TextBlockMarked())
					DeleteSelectedText();
				else 	{
	            	DeleteCharacter();
					PaintCurrentLine();
				}
            	return;
			case INS:
				insertmode = insertmode ? False : True;
				ResetCursor();
				break;
			case CTRL_X:
				Cut();
				return;
			case CTRL_C:
				Copy();
				return;
			case CTRL_V:
				Paste();
				return;
        	default:
            	if (!isprint(key))
                	break;
            	// --- all printable keys get processed by editbox
				InsertCharacter(key);
				BuildTextPointers();
				PaintCurrentLine();
            	return;
    	}
	}
	TextBox::Keyboard(key);
}

// -------- move the editbox
void EditBox::Move(int x, int y)
{
	TextBox::Move(x, y);
	ResetCursor();
}

// --------- clear the text from the editbox
void EditBox::ClearText()
{
	TextBox::ClearText();
	OpenWindow();
	ResetCursor();
}

void EditBox::Home()
{
	column = 0;
	wleft = 0;
	Paint();
	ResetCursor();
}

void EditBox::End()
{
	int ch;
	while ((ch = CurrentChar()) != '\0' && ch != '\n')
		column++;
	if (column - wleft >= ClientWidth())	{
		wleft = column - ClientWidth() + 1;
		Paint();
	}
	ResetCursor();
}

// ------- page right one screenwidth
Bool EditBox::PageRight()
{
	if (TextBox::PageRight())	{
		column = wleft;
		ResetCursor();
		return True;
	}
	return False;
}

// ------- page left one screenwidth
Bool EditBox::PageLeft()
{
	if (TextBox::PageLeft())	{
		column = wleft;
		ResetCursor();
		return True;
	}
	return False;
}

void EditBox::NextWord()
{
	while (!isWhite(CurrentChar()) && CurrentChar())
		Forward();
	while (isWhite(CurrentChar()))
		Forward();
}

void EditBox::PrevWord()
{
	Backward();
	while (isWhite(CurrentChar()) && !AtBufferStart())
		Backward();
	while (!isWhite(CurrentChar()) && !AtBufferStart())
		Backward();
	if (isWhite(CurrentChar()))
		Forward();
}

void EditBox::Forward()
{
	if (CurrentChar())	{
		column++;
		if (column-wleft == ClientWidth())
			ScrollLeft();
		ResetCursor();
	}
}

void EditBox::Backward()
{
	if (column)	{
		if (column == wleft)
			ScrollRight();
		--column;
		ResetCursor();
	}
}

void EditBox::InsertCharacter(int key)
{
	DeleteSelectedText();
	unsigned col = CurrentCharPosition();
	if (textlength && text->Strlen() == textlength+1)	{
		desktop.speaker().Beep();
		return;
	}
	if (insertmode || CurrentChar() == '\0')	{
		// ---- shift the text to make room for new character
		String ls, rs;
		if (col)
			ls = text->left(col);
		int rt = text->Strlen()-col;
		if (rt > 0)
			rs = text->right(rt);
		*text = ls + " " + rs;
	}
	(*text)[col] = (char) key;
	if (key == '\n')
		BuildTextPointers();
	Forward();
	changed = True;
}

void EditBox::DeleteCharacter()
{
	if (CurrentChar())	{
		String ls, rs;
		unsigned col = CurrentCharPosition();
		if (col)
			ls = text->left(col);
		int rt = text->Strlen()-col-1;
		if (rt > 0)
			rs = text->right(rt);
		*text = ls + rs;
		BuildTextPointers();
		changed = True;
	}
}
// ---- position the cursor
void EditBox::SetCursor(int x, int y)
{
	column = x;
	desktop.cursor().SetPosition(
		x+ClientLeft()-wleft, y+ClientTop()-wtop);
}

void EditBox::LeftButton(int mx, int my)
{
	if (ClientRect().Inside(mx, my))	{
		column = max(0, min(text->Strlen()-1, mx-ClientLeft()+wleft));
		ResetCursor();
	}
	TextBox::LeftButton(mx, my);
}

void EditBox::SetCursorSize()
{
	if (insertmode)
	    desktop.cursor().BoxCursor();
	else 
	    desktop.cursor().NormalCursor();
}

Bool EditBox::ResetCursor()
{
	SetCursor(column, GetRow());
	if (isVisible() && this == desktop.InFocus())	{
		if (column >= wleft && column < wleft+ClientWidth())	{
			SetCursorSize();
		    desktop.cursor().Show();
			return True;
		}
	    desktop.cursor().Hide();
	}
	return False;
}

void EditBox::ExtendBlock(int x, int y)
{
	column = x;
	TextBox::ExtendBlock(x, y);
}

void EditBox::DeleteSelectedText()
{
	if (TextBlockMarked())	{
		if (BlkEndLine < BlkBegLine)
			column = BlkEndCol;
		else if (BlkEndLine > BlkBegLine)
			column = BlkBegCol;
		else
	        column = min(BlkBegCol, BlkEndCol);
        if (column < wleft || column >= wleft + ClientWidth())
            wleft = column;
		TextBox::DeleteSelectedText();
		changed = True;
	}
}

void EditBox::InsertText(const String& txt)
{
	int ln1 = CurrentPosition();
	int ln2 = TextLength() - ln1;
	*text = text->left(ln1) + txt + text->right(ln2);
	BuildTextPointers();
	Paint();
}

void EditBox::Cut()
{
	if (TextBlockMarked())	{
		desktop.clipboard().SetText(GetSelectedText());
		DeleteSelectedText();
	}
}

void EditBox::Copy()
{
	if (TextBlockMarked())
		desktop.clipboard().SetText(GetSelectedText());
}

void EditBox::Paste()
{
	InsertText(desktop.clipboard().GetText());
}


