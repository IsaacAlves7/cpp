// ----- editor.cpp

#include "editor.h"
#include "desktop.h"

// ----------- common constructor code
void Editor::OpenWindow()
{
    windowtype = EditorWindow;
    row = 0;
    tabs = 4;
    insertmode = desktop.keyboard().InsertMode();
    wordwrapmode = True;
	SetSingleBorder();
}
// ---- keep the cursor out of tabbed space
void Editor::AdjustCursorTabs(int key)
{
    while (CurrentChar() == Ptab)
        key == FWD ? column++ : --column;
    ResetCursor();
}
// -------- process keystrokes
void Editor::Keyboard(int key)
{
    int svwtop = wtop;
    int svwleft = wleft;
    switch (key)    {
        case '\t':
            InsertTab();
            BuildTextPointers();
            PaintCurrentLine();
            ResetCursor();
            break;
        case ALT_P:
			ClearTextBlock();
            FormParagraph();
            break;
        case UP:
            Upward();
            TestMarking();
            break;
        case DN:
            Downward();
            TestMarking();
            break;
        case CTRL_HOME:
            BeginDocument();
            TestMarking();
            break;
        case CTRL_END:
            EndDocument();
            TestMarking();
            break;
        case '\r':
            InsertCharacter('\n');
            BuildTextPointers();
            ResetCursor();
            Paint();
            break;
        case DEL:
        case RUBOUT:
            EditBox::Keyboard(key);
			WordWrap();
            break;
        default:
            EditBox::Keyboard(key);
            break;
    }
    if (svwtop != wtop || svwleft != wleft)
        Paint();
}
// --- move the cursor forward one character
void Editor::Forward()
{
    if (CurrentChar())    {
        if (CurrentChar() == '\n')    {
            Home();
            Downward();
        }
        else
            EditBox::Forward();
        AdjustCursorTabs(FWD);
    }
}
// --- move the cursor back one character
void Editor::Backward()
{
    if (column)
        EditBox::Backward();
    else if (row)    {
        Upward();
        End();
    }
    AdjustCursorTabs();
}
// ---- if cursor moves out of the window, scroll
void Editor::ScrollCursor()
{
    if (column < wleft || column >= wleft + ClientWidth())    {
        wleft = column;
        Paint();
    }
}
// --- move the cursor up one line
void Editor::Upward()
{
    if (row)    {
        if (row == wtop)
            ScrollDown();
        --row;
        AdjustCursorTabs();
        ScrollCursor();
    }
}
// --- move the cursor down one line
void Editor::Downward()
{
    if (row < wlines)    {
        if (row == wtop + ClientHeight() - 1)
            ScrollUp();
        row++;
        AdjustCursorTabs();
        ScrollCursor();
    }
}
// --- move the cursor to the beginning of the document
void Editor::BeginDocument()
{
    row = 0;
    wtop = 0;
    EditBox::Home();
    AdjustCursorTabs();
}
// --- move the cursor to the end of the document
void Editor::EndDocument()
{
    TextBox::End();
    row = wlines-1;
    End();
    AdjustCursorTabs();
}
// --- keep cursor in the window, in text and out of empty space
Bool Editor::ResetCursor()
{
    KeepInText(column, row);
    if (EditBox::ResetCursor())    {
        if (!(row >= wtop && row < wtop+ClientHeight()))    {
            desktop.cursor().Hide();
            return False;
        }
    }
    return True;
}
// ------- page up one screenfull
Bool Editor::PageUp()
{
    if (wlines)    {
        row -= ClientHeight();
        if (row < 0)
            row = 0;
        EditBox::PageUp();
        AdjustCursorTabs();
        return True;
    }
    return False;
}
// ------- page down one screenfull
Bool Editor::PageDown()
{
    if (wlines)    {
        row += ClientHeight();
        if (row >= wlines)
            row = wlines-1;
        EditBox::PageDown();
        AdjustCursorTabs();
        return True;
    }
    return False;
}
// --- insert a tab into the edit buffer
void Editor::InsertTab()
{
    ClearVisible();
    if (insertmode)    {
        EditBox::InsertCharacter('\t');
        while ((column % tabs) != 0)
            EditBox::InsertCharacter(Ptab);
    }
    else
        do
            Forward();
        while ((column % tabs) != 0);
    SetVisible();
}
// --- When inserting char, adjust next following tab, same line
void Editor::AdjustTabInsert()
{
    ClearVisible();
    // ---- test if there is a tab beyond this character
    int savecol = column;
    while (CurrentChar() && CurrentChar() != '\n')    {
        if (CurrentChar() == '\t')    {
            column++;
            if (CurrentChar() == Ptab)
                EditBox::DeleteCharacter();
            else
                for (int i = 0; i < tabs-1; i++)
                    EditBox::InsertCharacter(Ptab);
            break;
        }
        column++;
    }
    column = savecol;
    SetVisible();
}
// --- test for wrappable word and wrap it
void Editor::WordWrap()
{
    // --- test for word wrap
    int len = LineLength(row);
    int wd = ClientWidth()-1;
    if (len >= wd)    {
        const char *cp = TextLine(row);
        char ch = *(cp + wd);
        // --- test words beyond right margin
        if (len > wd || (ch && ch != ' ' && ch != '\n'))    {
            // --- test typing in last word in window's line
            const char *cw = cp + wd;
            cp += column;
            while (cw > cp)    {
                if (*cw == ' ')
                    break;
                --cw;
            }
            int newcol = 0;
            if (cw <= cp)    {
                // --- user was typing last word on line
                // --- find beginning of the word
                const char *cp1 = TextLine(row);
                const char *cw1 = cw;
                while (*cw1 != ' ' && cw1 > cp1)
                    --cw1, newcol++;
                wleft = 0;
            }
            FormParagraph();
            if (cw <= cp)    {
                // --- user was typing last word on line
                column = newcol;
                if (cw == cp)
                    --column;
                row++;
                if (row - wtop >= ClientHeight())
                    ScrollUp();
                ResetCursor();
            }
        }
    }
	else if (row < wlines-2)	{
		// ---- room at the end of the line and at least one more line
		int room = wd-len;
		// ---- measure the first word on the next line
        const char *cp = TextLine(row+1);
		if (*cp != '\t')	{
			// --- bypass and count white space
			int wh = 0;
			while (wh < tabs && *cp == ' ' && (*cp & 0x7f) == '\t')
				cp++, wh++;
			if (wh < tabs && *cp != '\n' && *cp)	{
				// --- found the first word, measure it
				int wdct = 0;
				while (*cp != ' ' && *cp != '\n' && (*cp & 0x7f) != '\t')
					wdct++, cp++;
				// --- if room on this line for 1st word on next...
				if (wdct < room)
					FormParagraph();
			}
		}
	}
}
// --- insert a character at the current cursor position
void Editor::InsertCharacter(int key)
{
    if (insertmode)    {
        if (key != '\n')
            AdjustTabInsert();
    }
    else if (CurrentChar() == '\t')    {
        // --- overtyping a tab
        ClearVisible();
        column++;
        while (CurrentChar() == Ptab)
            EditBox::DeleteCharacter();
        --column;
    }
    ClearVisible();
    EditBox::InsertCharacter(key);
    SetVisible();
    ResetCursor();
    if (wordwrapmode)
        WordWrap();
}
// --- When deleting char, adjust next following tab, same line
void Editor::AdjustTabDelete()
{
    ClearVisible();
    // ---- test if there is a tab beyond this character
    int savecol = column;
    while (CurrentChar() && CurrentChar() != '\n')    {
        if (CurrentChar() == '\t')    {
            column++;
            // --- count pseudo tabs
            int pct = 0;
            while (CurrentChar() == Ptab)
                pct++, column++;
            if (pct == tabs-1)    {
                column -= tabs-1;
                for (int i = 0; i < tabs-1; i++)
                    EditBox::DeleteCharacter();
            }
            else
                EditBox::InsertCharacter(Ptab);
            break;
        }
        column++;
    }
    column = savecol;
    SetVisible();
}
// --- delete the character at the current cursor position
void Editor::DeleteCharacter()
{
    if (CurrentChar() == '\0')
        return;
    if (insertmode)
        AdjustTabDelete();
    if (CurrentChar() == '\t')    {
        // --- deleting a tab
        EditBox::DeleteCharacter();
        while (CurrentChar() == Ptab)
            EditBox::DeleteCharacter();
        return;
    }
    const char *cp = TextLine(row);
    const char *cw = cp + column;
    const Bool delnewline = (Bool) (*cw == '\n');
    const Bool reform = (Bool) (delnewline && *(cw+1) != '\n' && column);
    const Bool lastnewline =
                        (Bool) (delnewline && *(cw+1) == '\0');
    int newcol = 0;
    if (reform && !lastnewline)
        // --- user is deleting /n, find beginning of last word
        while (cw > cp && *--cw != ' ')
            newcol++;
    EditBox::DeleteCharacter();
    if (!lastnewline)	{
    	if (delnewline && !reform)
        	// --- user deleted a blank line
        	Paint();
    	else if (wordwrapmode && reform)    {
        	// --- user deleted /n
        	wleft = 0;
        	FormParagraph();
        	if (CurrentChar() == '\n')    {
            	// ---- joined the last word with next line's
            	//      first word and then wrapped the result
            	column = newcol;
            	row++;
        	}
    	}
	}
}
// --- form a paragraph from the current cursor position
//    through one line before the next blank line or end of text
void Editor::FormParagraph()
{
    int BegCol, FirstLine;
    const char *blkBegLine, *blkEndLine, *blkBeg;

    // ---- forming paragraph from cursor position
    FirstLine = wtop + row;
    blkBegLine = blkEndLine = TextLine(row);
    if ((BegCol = column) >= ClientWidth())
        BegCol = 0;
    // ---- locate the end of the paragraph
    while (*blkEndLine)    {
        Bool newpara = True;
		int x = 0;
        // --- leading tab or spaces mark end of paragraph
        while (*(blkEndLine+x) && *(blkEndLine+x) != '\n')    {
			// --- test \t in left margin
			if (*blkEndLine != '\t')
				// --- test spaces in left margin
	            if (x < tabs && *(blkEndLine+x) != ' ')
    	            newpara = False;
			x++;
        }
        if (newpara)
            break;
        blkEndLine += x;
        if (*blkEndLine)
			blkEndLine++;
    }
    --blkEndLine;
    if (blkEndLine == blkBegLine)    {
        SetVisible();
        Downward();
        return;
    }
    // --- change newlines, tabs, and tab expansions to spaces
    blkBeg = blkBegLine;
    while (blkBeg < blkEndLine)    {
        if (*blkBeg == '\n' || ((*blkBeg) & 0x7f) == '\t')    {
            int off = blkBeg - (const char *)*text;
            (*text)[off] = ' ';
        }
        blkBeg++;
    }
    // ---- insert newlines at new margin boundaries
    blkBeg = blkBegLine;
    while (blkBegLine < blkEndLine)    {
        blkBegLine++;
        if ((int)(blkBegLine - blkBeg) == ClientWidth()-1)    {
            while (*blkBegLine != ' ' && blkBegLine > blkBeg)
                --blkBegLine;
            if (*blkBegLine != ' ')    {
                blkBegLine = strchr(blkBegLine, ' ');
				while (*blkBegLine && *blkBegLine != ' ')
					blkBegLine++;
                if (*blkBegLine == '\0' ||
                        blkBegLine >= blkEndLine)
                    break;
            }
            int off = blkBegLine - (const char *)*text;
            (*text)[off] = '\n';
            blkBeg = blkBegLine+1;
        }
    }
    BuildTextPointers();
    changed = True;
    // --- put cursor back at beginning
    column = BegCol;
    if (FirstLine < wtop)
        wtop = FirstLine;
    row = FirstLine - wtop;
    SetVisible();
    Paint();
    ResetCursor();
}
// --------- add a line of text to the editor textbox
void Editor::AddText(const String& txt)
{
    // --- compute the buffer size based on tabs in the text
    const char *tp = txt;
    int x = 0;
    int sz = 0;
    while (*tp)    {
        if (*tp == '\t')    {
            // --- tab, adjust the buffer length
            int sps = Tabs() - (x % Tabs());
            sz += sps;
            x += sps;
        }
        else    {
            // --- not a tab, count the character
            sz++;
            x++;
        }
        if (*tp == '\n')
            x = 0;    // newline, reset x
        tp++;
    }
    // --- allocate a buffer
    char *ep = new char[sz];
    // --- detab the input file
    tp = txt;
    char *ttp = ep;
    x = 0;
    while (*tp)    {
        // --- put the character (\t, too) into the buffer
        *ttp++ = *tp;
        x++;
        // --- expand tab into \t and expansions (\t + 0x80)
        if (*tp == '\t')
            while ((x % Tabs()) != 0)
                *ttp++ = Ptab, x++;
        else if (*tp == '\n')
            x = 0;
        tp++;
    }
    *ttp = '\0';
    // ---- add the text to the editor window
    EditBox::AddText(String(ep));
}
// ------- retrieve editor text collapsing tabs
const String Editor::GetText()
{
	// --- build a char copy, collapsing tabs
    char *tx = new char[text->Strlen()+1];
    const char *tp = (const char *) *text;
    char *nt = tx;
    while (*tp)    {
        if (*(const unsigned char *)tp != Ptab)
            *tx++ = *tp;
        tp++;
    }
    *tx = '\0';
	// --- return the collapsed version of the text
    String temp(nt);
    return nt;
}
// --- write a string to the editor window
void Editor::WriteString(const String& ln,int x,int y,int fg,int bg)
{
    String *nln = new String(ln.Strlen(), 0);
    int ch;
    for (int i = 0; i < ln.Strlen(); i++)    {
        ch = ln[i];
        (*nln)[i] = (ch & 0x7f) == '\t' ? ' ' : ch;
    }
    EditBox::WriteString(*nln, x, y, fg, bg);
	delete nln;
}
// ---- left mouse button pressed
void Editor::LeftButton(int mx, int my)
{
    if (ClientRect().Inside(mx, my))    {
        column = mx-ClientLeft()+wleft;
        row = my-ClientTop()+wtop;
        ResetCursor();
    }
    TextBox::LeftButton(mx, my);
}
// --------- clear the text from the editor window
void Editor::ClearText()
{
    row = 0;
    EditBox::ClearText();
}
// ------ extend the marked block
void Editor::ExtendBlock(int x, int y)
{
    row = y;
    EditBox::ExtendBlock(x, y);
}
// ---- delete a marked block
void Editor::DeleteSelectedText()
{
    if (TextBlockMarked())    {
        row = min(BlkBegLine, BlkEndLine);
        if (row < wtop || row >= wtop + ClientHeight())
            wtop = row;
        EditBox::DeleteSelectedText();
        FormParagraph();
    }
}
// ---- insert a string into the text
void Editor::InsertText(const String& txt)
{
    EditBox::InsertText(txt);
    FormParagraph();
}
// ---- set tab width
void Editor::SetTabs(int t)
{
    if (t && tabs != t && t <= MaxTab)    {
        tabs = t;
        if (text != 0)    {
            // ------- retab the text
            String *ln = new String;
            for (int lno = 0; lno < wlines; lno++)    {
                // --- retrieve a line at a time
                ExtractTextLine(*ln, lno);
                int len = ln->Strlen();
                String *newln = new String(len * t, 0);
                // --- copy string, collapsing old tabs
                //     and expanding new tabs
                for (int x2 = 0, x1 = 0; x2 < len; x2++)    {
                    unsigned int ch = (*ln)[x2] & 0xff;
                    if (ch == Ptab)
                        // --- collapse old tab expansion
                        continue;
                    // --- copy text
                    (*newln)[x1++] = ch;
                    if (ch == '\t')
                        // --- expand new tabs
                        while ((x1 % t) != 0)
                            (*newln)[x1++] = Ptab;
                }
                (*newln)[x1] = '\n';
				newln->ChangeLength(x1+1);

                // --- compute left segment length
                unsigned seg1 = (unsigned)
                    (TextLine(lno) - (const char *) *text);
                // --- compute right segment length
                unsigned seg2 = 0;
                if (lno < wlines-1)
                    seg2 = (unsigned) text->Strlen() - 
                    (TextLine(lno+1) - (const char *) *text);
                // --- rebuild the text from the three parts
                String *lft = new String(text->left(seg1));
                String *rht = new String(text->right(seg2));
                *text = *lft + *newln + *rht;

				delete rht;
				delete lft;
				delete newln;
                BuildTextPointers();
            }
			delete ln;
            Paint();
        }
    }
}
// ---- position the cursor
void Editor::SetCursor(int x, int y)
{
	row = y;
	EditBox::SetCursor(x, y);
}


