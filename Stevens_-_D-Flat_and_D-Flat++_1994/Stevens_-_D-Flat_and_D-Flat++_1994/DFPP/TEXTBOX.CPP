// ------------- textbox.cpp

#include "desktop.h"
#include "textbox.h"
#include "scrolbar.h"

static Color col = {
	BLACK,				// fg
	LIGHTGRAY,			// bg
	LIGHTGRAY,			// selected fg
	BLACK,				// selected bg
	LIGHTGRAY,			// frame fg
	BLACK,				// frame bg
	LIGHTGRAY,			// highlighted fg
	BLUE				// highlighted bg
};

// ----------- common constructor code
void TextBox::OpenWindow()
{
    windowtype = TextboxWindow;
    text = 0;
    hscrollbar = vscrollbar = 0;
    textpointers = 0;
	SetColor(col);
    ClearText();
	textlength = 0;
}

void TextBox::CloseWindow()
{
    ClearText();
    delete hscrollbar;
    delete vscrollbar;
    Control::CloseWindow();
}


// ------ show the textbox
void TextBox::Show()
{
    if ((Attribute() & HSCROLLBAR) && hscrollbar == 0)    {
        hscrollbar = new ScrollBar(HORIZONTAL, this);
        hscrollbar->SetAttribute(FRAMEWND);
    }
    if ((Attribute() & VSCROLLBAR) && vscrollbar == 0)    {
        vscrollbar = new ScrollBar(VERTICAL, this);
        vscrollbar->SetAttribute(FRAMEWND);
    }
    Control::Show();
}

// ------------ build the text line pointers
void TextBox::BuildTextPointers()
{
    textwidth = wlines = 0;
    // ---- count the lines of text
	const char *cp = *text;
    while (*cp)    {
        wlines++;
        while (*cp && *cp != '\n')
            cp++;
        if (*cp)
            cp++;
    }
    // ----- build the pointer array
    delete textpointers;
    textpointers = new unsigned int[wlines+1];
	*textpointers = 0;
    unsigned int off = 0;
    cp = *text;
    wlines = 0;
    while (*cp)    {
        *(textpointers + wlines++) = off;
        const char *cp1 = cp;
        while (*cp && *cp != '\n')
            cp++;
        textwidth = max((unsigned int) textwidth, (unsigned int) (cp - cp1));
        if (*cp)
            cp++;
        off = (unsigned int) (cp - *text);
    }
}

// --------- add a line of text to the textbox
void TextBox::AddText(const String& txt)
{
	int off = 0;
    if (text == 0)
		// --- first line being added
        text = new String(txt);
	else	{
		// --- line added to existing lines
		off = text->Strlen();
	    *text += txt;
	}
	// ---- assure that text ends with a newline
	if (txt[txt.Strlen()-1] != '\n')
		*text += "\n";

	// --- add an entry to the textpointers offset array
    unsigned *tp = new unsigned int[wlines+1];
	if (textpointers)	{
		for (int i = 0; i < wlines; i++)
			*(tp+i) = *(textpointers+i);
		delete textpointers;
	}
	textpointers = tp;
    *(textpointers + wlines++) = off;

	// --- assure that textbox does not exceed its permitted length
	if (textlength && text->Strlen() > textlength)	{
		text->ChangeLength(textlength);
		BuildTextPointers();
	}
}

// --------- set the textbox's text buffer to new text
void TextBox::SetText(const String& txt)
{
    ClearText();
	AddText(txt);
	BuildTextPointers();
}

// ------ set the length of the text buffer
void TextBox::SetTextLength(unsigned int len)
{
    if (len > 0 && text != 0)
        text->ChangeLength(len);
	textlength = len;
}

// --------- clear the text from the textbox
void TextBox::ClearText()
{
    delete text;
	text = 0;
    wlines = 0;
    textwidth = 0;
    wtop = wleft = 0;
    ClearTextBlock();
    delete textpointers;
	textpointers = 0;
}

void TextBox::ClearTextBlock()
{
	BlkBegLine=BlkEndLine=BlkBegCol=BlkEndCol=0;
	ancx = ancy = 0;
	buttondown = False;
	textmarking = NotMarking;
}

// ----- mark a block of selected text
void TextBox::SelectText(int lf, int tp, int rt, int bt)
{
	BlkBegCol = lf;
	BlkBegLine = tp;
	BlkEndCol = rt;
	BlkEndLine = bt;
}

// ------- extract a text line
void TextBox::ExtractTextLine(String& ln, int lno)
{
    const char *lp = TextLine(lno);
    int offset = lp - (const char *) *text;
    for (int len = 0; *(lp+len) && *(lp+len) != '\n'; len++)
        ;
    ln = text->mid(len, offset);
}

// ---- display a line with a shortcut key character
void TextBox::WriteShortcutLine(int lno, int fg, int bg)
{
    if (!isVisible())
        return;
    String *sc = new String;
	ExtractTextLine(*sc, lno);
    int x = sc->Strlen();
    int y = lno-wtop;
    x -= DisplayShortcutField(*sc, 0, y, fg, bg);
    // --------- pad the line
    int wd = ClientWidth() - x;
    if (wd > 0)
        WriteClientString(String(wd, ' '), x, y, fg, bg);
	delete sc;
}

// ---- display a shortcut field character
int TextBox::DisplayShortcutField(const String& sc, int x, int y,
                                                int fg, int bg)
{
    int scs = 0;
    if (isVisible())	{
    	int off = sc.FindChar(SHORTCUTCHAR);
    	if (off != -1)    {
        	scs++;
        	if (off != 0)    {
            	String *ls = new String(sc.left(off));
            	WriteClientString(*ls, x, y, fg, bg);
				delete ls;
        	}
        	WriteClientChar(sc[off+1], x+off, y, shortcutfg, bg);
        	int len = sc.Strlen()-off-2;
        	if (len > 0)    {
            	String *rs = new String(sc.right(len));
            	scs += DisplayShortcutField(*rs, x+off+1, y, fg, bg);
				delete rs;
        	}
    	}
    	else
        	WriteClientString(sc, x, y, fg, bg);
	}
    return scs;
}

// ------- write a text line to the textbox
void TextBox::WriteTextLine(int lno)
{
    if (!isVisible() || lno < wtop || lno >= wtop + ClientHeight())
        return;
    if (TextBlockMarked())    {
		int bbc = BlkBegCol;
		int bec = BlkEndCol;
		int bbl = BlkBegLine;
		int bel = BlkEndLine;

        // ----- put lowest marker first
        if (bbl > bel)    {
            swap(bbl, bel);
            swap(bbc, bec);
        }
        if (bbl == bel && bbc > bec)
            swap(bbc, bec);

        if (lno >= bbl && lno <= bel)    {
            // ------ the block includes this line
    		String *ln = new String;
			ExtractTextLine(*ln, lno);
			int wd = ClientWidth();
			int len = ln->Strlen();

			// --- compute length of 1st unmarked segment
			int u1len = 0;
			if (lno == bbl)
				u1len = bbc;

			// --- compute length of 2nd marked segment
			int mklen = len-u1len;
			if (lno == bel)
				mklen = bec-u1len;

			// --- compute length of 3rd unmarked segment
			int u3len = 0;
			if (lno == bel)
				u3len = len - (u1len + mklen);

			// --- now build the three strings
			String *fst = new String;
			String *mrk = new String;
			String *lst = new String;
			if (u1len > 0)
				*fst = ln->left(u1len);
			if (mklen > 0)
				*mrk = ln->mid(mklen, u1len);
			if (u3len > 0)
				*lst = ln->right(u3len);

			// --- truncate the strings according to window panning
			if (u1len > 0 && wleft > 0)	{
				u1len -= wleft;
				if (u1len > 0)
					*fst = fst->right(u1len);
				else
					u1len = 0;
			}

			if (mklen > 0 && bec >= wleft && bbc-wleft < wd)	{
				int bc = 0;
				if (bbl == lno)
					bc = bbc;
				if (bc < wleft)	{
					mklen -= wleft - bc;
					*mrk = mrk->right(mklen);
				}
				if (u1len+mklen > wd)	{
					mklen = wd-u1len;
					if (mklen > 0)
						*mrk = mrk->left(mklen);
				}
			}
			else
				mklen = 0;

			if (bbc-wleft >= wd)
				u3len = 0;
			if (u3len > 0)	{
				u3len = min(u3len, wd-(u1len+mklen));
				if (u3len > 0)
					*lst = lst->left(u3len);
			}

			u1len = max(0, u1len);
			mklen = max(0, mklen);
			u3len = max(0, u3len);

			// --- display the three strings
			if (u1len)
				WriteClientString(*fst, 0, lno-wtop, ClientFG(), ClientBG());
			if (mklen)
				WriteClientString(*mrk, u1len, lno-wtop, SelectedFG(), SelectedBG());
			int pd = mklen - mrk->Strlen();
			if (pd > 0)	{
				String *pad = new String(pd, ' ');
				int poff = u1len + mrk->Strlen();
				WriteClientString(*pad, poff, lno-wtop, SelectedFG(), SelectedBG());
				delete pad;
			}
			if (u3len)
				WriteClientString(*lst, u1len+mklen,	lno-wtop, ClientFG(), ClientBG());

			//  ---- pad the string display
			PadString(lno, u1len + mklen + u3len);

			delete lst;
			delete mrk;
			delete fst;
			delete ln;
			return;
		}	
	}
	// --- none of line is marked
	WriteTextLine(lno, ClientFG(), ClientBG());
}

// ---------- pad a display string
void TextBox::PadString(int lno, int datlen)
{
	int padlen = ClientWidth() - datlen;
	if (padlen > 0)	{
		String *pad = new String(padlen, ' ');
		WriteClientString(*pad, datlen, lno-wtop, ClientFG(), ClientBG());
		delete pad;
	}
}

// ------- write a text line to the textbox
void TextBox::WriteTextLine(int lno, int fg, int bg)
{
    if (!isVisible() || lno < wtop || lno >= wtop + ClientHeight())
        return;
    int wd = ClientWidth();
    String *ln = new String;
	ExtractTextLine(*ln, lno);
	if (wleft)
		*ln = ln->mid(wd, wleft);
    // ----- display the line
    WriteClientString(*ln, 0, lno-wtop, fg, bg);
	// ----- display the line's padding
	PadString(lno, ln->Strlen());
}

// ---------- paint the textbox
void TextBox::Paint()
{
    if (isVisible())    {
    	if (text == 0 || wlines == 0)
        	Control::Paint();
    	else    {
        	int ht = ClientHeight();
        	int wd = ClientWidth();
        	for (int i = 0; i < min(wlines-wtop,ht); i++)
            	WriteTextLine(wtop+i);
        	// ---- pad the bottom lines in the window
        	String *line = new String(wd, ' ');
        	while (i < ht)
            	WriteClientString(*line, 0, i++, ClientFG(), ClientBG());
			delete line;
        	if (resetscrollbox)
            	SetScrollBoxes();
        	resetscrollbox = False;
    	}
	}
}

// ------ process a textbox keystroke
void TextBox::Keyboard(int key)
{
    switch (key)    {
        case UP:
            if (ClientTop() == ClientBottom())
                break;
            ScrollDown();
            return;
        case DN:
            if (ClientTop() == ClientBottom())
                break;
            ScrollUp();
            return;
        case FWD:
            ScrollLeft();
            return;
        case BS:
            ScrollRight();
            return;
        case PGUP:
            PageUp();
			TestMarking();
            return;
        case PGDN:
            PageDown();
			TestMarking();
            return;
        case CTRL_PGUP:
            PageLeft();
			TestMarking();
            return;
        case CTRL_PGDN:
            PageRight();
			TestMarking();
            return;
        case HOME:
            Home();
            return;
        case END:
            End();
            return;
        default:
            break;
    }
    Control::Keyboard(key);
}

// ------- scroll up one line
Bool TextBox::ScrollUp()
{
    if (wtop < wlines-1)    {
        desktop.screen().Scroll(ClientRect(), 1, ClientFG(), ClientBG());
        wtop++;
        int ln = wtop+ClientHeight()-1;
        if (ln < wlines)
            WriteTextLine(ln);
        SetScrollBoxes();
		return True;
    }
	return False;
}

// ------- scroll down one line
Bool TextBox::ScrollDown()
{
    if (wtop)    {
        desktop.screen().Scroll(ClientRect(), 0, ClientFG(), ClientBG());
        --wtop;
        WriteTextLine(wtop);
        SetScrollBoxes();
		return True;
    }
	return False;
}

// ------- scroll left one character
Bool TextBox::ScrollLeft()
{
    if (wleft < textwidth-ClientWidth())	{
        wleft++;
	    Paint();
		return True;
	}
	return False;
}

// ------- scroll right one character
Bool TextBox::ScrollRight()
{
    if (wleft > 0)	{
        --wleft;
	    Paint();
		return True;
	}
	return False;
}

// ------- page up one screenfull
Bool TextBox::PageUp()
{
    if (wtop)    {
        wtop -= ClientHeight();
        if (wtop < 0)
            wtop = 0;
        resetscrollbox = True;
        Paint();
		return True;
    }
	return False;
}

// ------- page down one screenfull
Bool TextBox::PageDown()
{
    if (wtop < wlines-1)    {
        wtop += ClientHeight();
        if (wlines < wtop)
            wtop = wlines-1;
        resetscrollbox = True;
        Paint();
		return True;
    }
	return False;
}

// ------- page right one screenwidth
Bool TextBox::PageRight()
{
    if (wleft < textwidth-1)    {
        wleft += ClientWidth();
        if (wleft >= textwidth)
            wleft = textwidth-ClientWidth();
		if (wleft < 0)
			wleft = 0;
        resetscrollbox = True;
        Paint();
		return True;
    }
	return False;
}

// ------- page left one screenwidth
Bool TextBox::PageLeft()
{
    if (wleft)    {
        wleft -= ClientWidth();
        if (wleft < 0)
            wleft = 0;
        resetscrollbox = True;
        Paint();
		return True;
    }
	return False;
}

// ----- move to the first line of the textbox
void TextBox::Home()
{
    wtop = 0;
    Paint();
}

// ----- move to the last line of the textbox
void TextBox::End()
{
    wtop = wlines-ClientHeight();
    if (wtop < 0)
        wtop = 0;
    Paint();
}

// ----- position the scroll boxes
void TextBox::SetScrollBoxes()
{
    if (vscrollbar != 0)
        vscrollbar->TextPosition(wlines ? (wtop*100)/wlines : 0);
    if (hscrollbar != 0)
        hscrollbar->TextPosition(textwidth ? (wleft*100)/textwidth : 0);
}

// ---- compute the horizontal page position
void TextBox::HorizontalPagePosition(int pct)
{
    wleft = (textwidth * pct) / 100;
    Paint();
}

// ---- compute the vertical page position
void TextBox::VerticalPagePosition(int pct)
{
    wtop = (wlines * pct) / 100;
    Paint();
}

void TextBox::SetScrollBars()
{
	if (LineCount() > ClientHeight())
		SetAttribute(VSCROLLBAR);
	else	{
		ClearAttribute(VSCROLLBAR);
		delete vscrollbar;
		vscrollbar = 0;
	}
	if (TextWidth() > ClientWidth())
		SetAttribute(HSCROLLBAR);
	else	{
		ClearAttribute(HSCROLLBAR);
		delete hscrollbar;
		hscrollbar = 0;
	}
}

// ---- determine the length of a line of text
int TextBox::LineLength(int lno)
{
	const char *tp = TextLine(lno);
	int len = 0;
	while (*(tp + len) && *(tp + len) != '\n')
		len++;
	return len;
}

void TextBox::SetAnchor(int x, int y)
{
	BlkBegLine = y;
	BlkBegCol = x;
	BlkEndLine = y;
	BlkEndCol = x;
}

void TextBox::ExtendBlock(int x, int y)
{
	int oldline = BlkEndLine;
	BlkEndLine = max(0, min(wlines-1, y));
	BlkEndCol = max(0, min(LineLength(BlkEndLine), x));
	int newline = BlkEndLine;
	if (newline < oldline)
		swap(newline, oldline);
	while (oldline <= newline)
		WriteTextLine(oldline++);
}

void TextBox::KeepInText(int& col, int& row)
{
    row = min(row, wlines);
    col = min(col, LineLength(row));
}

void TextBox::LeftButton(int mx, int my)
{
	if (ClientRect().Inside(mx, my))	{
		if (textmarking == MouseMarking)	{
			// --- marking text and landed in border?
			if (mx == ClientLeft())	{
				if (ScrollRight())
					ExtendBlock(BlkEndCol-1, BlkEndLine);
			}
			else if (mx == ClientRight())	{
				if (ScrollLeft())
					ExtendBlock(BlkEndCol+1, BlkEndLine);
			}
			else if (my == ClientTop())	{
				if (ScrollDown())
					ExtendBlock(BlkEndCol, BlkEndLine-1);
			}
			else if (my == ClientBottom())	{
				if (ScrollUp())
					ExtendBlock(BlkEndCol, BlkEndLine+1);
			}
		}
		else	{
			if (TextBlockMarked())	{
				ClearTextBlock();
				Paint();
			}
			ancx = mx-ClientLeft()+wleft;
			ancy = my-ClientTop()+wtop;
			KeepInText(ancx, ancy);
			buttondown = True;
		}
	}
	else
		Control::LeftButton(mx, my);
}

void TextBox::MouseMoved(int mx, int my)
{
	if (ClientRect().Inside(mx, my))	{
		int msx = mx-ClientLeft()+wleft;
		int msy = my-ClientTop()+wtop;
		KeepInText(msx, msy);
		if (buttondown)	{
			buttondown = False;
			SetAnchor(ancx, ancy);
			textmarking = MouseMarking;
			desktop.mouse().SetTravel(ClientLeft(),
									  ClientRight(),
									  ClientTop(),
									  ClientBottom());
		}
		if (textmarking == MouseMarking)
			ExtendBlock(msx, msy);
	}
	Control::MouseMoved(mx, my);
}

void TextBox::ButtonReleased(int mx, int my)
{
	buttondown = False;
	if (textmarking == MouseMarking)	{
		textmarking = NotMarking;
		desktop.mouse().RestoreTravel();
	}
	Control::ButtonReleased(mx, my);
}

void TextBox::ShiftChanged(int sk)
{
	Bool shftd = (Bool) (sk & (LEFTSHIFT | RIGHTSHIFT));
	if (shftd && textmarking != KeyboardMarking)	{
		ancx = GetColumn();
		ancy = GetRow();
	}
}

void TextBox::TestMarking()
{
	Bool shftd = (Bool) (desktop.keyboard().GetShift() &
					(LEFTSHIFT | RIGHTSHIFT));
	if (TextBlockMarked() && !shftd)	{
		ClearTextBlock();
		Paint();
	}
	else if (textmarking == NotMarking && shftd)	{
		SetAnchor(ancx, ancy);
		textmarking = KeyboardMarking;
	}
	if (textmarking && shftd)
		ExtendBlock(GetColumn(), GetRow());
}

int TextBox::BlockBeginOffset() const
{
	int bb = (int) (TextLine(BlkBegLine) - *text) + BlkBegCol;
	int be = (int) (TextLine(BlkEndLine) - *text) + BlkEndCol;
	return min(bb, be);
}

int TextBox::BlockEndOffset() const
{
	int bb = (int) (TextLine(BlkBegLine) - *text) + BlkBegCol;
	int be = (int) (TextLine(BlkEndLine) - *text) + BlkEndCol;
	return max(bb, be);
}

String TextBox::GetSelectedText() const
{
	String stext;
	if (text != 0 && TextBlockMarked())
		stext = text->mid(BlockLength(), BlockBeginOffset());
	return stext;
}

void TextBox::DeleteSelectedText()
{
	if (text != 0 && TextBlockMarked())	{
		int ln1 = BlockBeginOffset();
		int ln2 = TextLength() - (BlockLength() + ln1);
		*text = text->left(ln1) + text->right(ln2);
		ClearTextBlock();
		BuildTextPointers();
		Paint();
	}
}


