// ------------- ted.cpp

#include <fstream.h>
#include <strstrea.h>
#include "ted.h"

String untitled("(untitled)");

extern unsigned _stklen = 8192;

int main(int argc, char *argv[])
{
    Ted ma;
	if (argc > 1)
		ma.OpenFile(argv[1]);
    ma.Execute();
    return 0;
}

// ------- editor colors
static Color ecol = {
	BLACK,				// fg
	LIGHTGRAY,			// bg
	LIGHTGRAY,			// selected fg
	BLACK,				// selected bg
	LIGHTGRAY,			// frame bg
	BLUE,				// frame fg
	LIGHTGRAY,			// highlighted fg
	BLUE				// highlighted bg
};

// ----- application colors
static Color acol = {
	LIGHTGRAY, BLUE,
	LIGHTGRAY, BLUE,
	LIGHTGRAY, BLUE,
	BLACK,	   CYAN		// ruler colors
};

// ---- construct application
Ted::Ted() :  menubar(TedMenu, this),
              toolbar(this),
              editor(ClientLeft(),
                     ClientTop(),
                     ClientHeight(),
                     ClientWidth(),
                     this),
              fname(untitled)
                 
{
	ruler = False;
	SetColor(acol);
    SetAttribute(SIZEABLE | MOVEABLE);
	editor.SetAttribute(VSCROLLBAR | BORDER);
	editor.SetColor(ecol);
    SetClearChar(' ');
    BuildTitle();
	if (!desktop.mouse().Installed())	{
		toolbar.Disable();
		ToolbarCmd.ClearToggle();
		CmToolbar();
	}
    Show();
	ruler = RulerCmd.isToggled();
	CmRuler();
	editor.SetFocus();
	desktop.SetHeapReserve(1024);
}

// ---- builds the title with the current document name
void Ted::BuildTitle()
{
    SetTitle(String("TED: ") + fname);
    Title();
}

// ---- File/New Menu Command
void Ted::CmNew()
{
    TestChanged();
    editor.ClearText();
	fname = untitled;
	BuildTitle();
    editor.Paint();
}

// ---- Open and load a specified file
void Ted::OpenFile(const String& fn)
{
    editor.ClearText();
    fname = fn;
    BuildTitle();
    editor.ClearChanged();

	NoticeBox *nb = new NoticeBox("Reading File");
    ifstream tfile(fname);
	if (tfile)	{
		// --- first get the file size
		tfile.seekg(0, ios::end);
		long sz = tfile.tellg();
		tfile.seekg(0, ios::beg);
		// --- buffer to read file into
    	char *ip = new char[sz+1];
		memset(ip, 0, sz+1);
		// --- read the file
		tfile.read(ip, sz);
		// --- put the text into the editor window
		String *buf = new String(ip);
		editor.SetText(*buf);
		// --- not needed any more
		delete buf;
		delete [] ip;
		delete nb;
    	editor.Paint();
	}
}

// ---- File/Open Menu Command
void Ted::CmOpen()
{
    TestChanged();
    FileOpen fo("*.txt");
    fo.Execute();
    if (fo.OKExit())
		OpenFile(fo.FileName());
}

// ---- File/Save Menu Command
void Ted::CmSave()
{
    if (fname == untitled)
        CmSaveAs();
    if (fname != untitled)    {
		NoticeBox("Saving File");
        editor.ClearChanged();
		int err;
        ofstream tfile;
		do	{
			tfile.open(fname);
    		err = desktop.TestCriticalError();
			if (err != 0)
				tfile.close();
		} while (err == 1);
		if (err == 0 && !tfile.fail())
	        tfile.write(editor.GetText(), editor.TextLength());
    }
}

// ---- File/Save As Menu Command
void Ted::CmSaveAs()
{
	String fs;
    if (fname != untitled)
		fs = fname;
    SaveAs sa("*.txt", fs);
    sa.Execute();
    if (sa.OKExit())    {
        fname = sa.FileName();
        BuildTitle();
        CmSave();
    }
}

// --- Prepare the edit menu
void Ted::PrepEdit()
{
	if (editor.TextBlockMarked())	{
		CopyCmd.Enable();
		CutCmd.Enable();
	}
	else	{
		CopyCmd.Disable();
		CutCmd.Disable();
	}
	if (desktop.clipboard().GetText().Strlen() > 0)
		PasteCmd.Enable();
	else 
		PasteCmd.Disable();
}

// ---- Edit/Cut Menu Command
void Ted::CmCut()
{
	editor.Cut();
}

// ---- Edit/Copy Menu Command
void Ted::CmCopy()
{
	editor.Copy();
}

// ---- Edit/Paste Menu Command
void Ted::CmPaste()
{
	editor.Paste();
}

// ---- Edit/Paragraph Menu Command
void Ted::CmPara()
{
	editor.FormParagraph();
}

// ---- View/Ruler Menu Command
void Ted::CmRuler()
{
	ruler = RulerCmd.isToggled();
	if (ruler)	{
		editor.Size(editor.Right(), editor.Bottom()-1);
		editor.Move(editor.Left(), editor.Top()+1);
		Paint();
	}
	else	{
		editor.Move(editor.Left(), editor.Top()-1);
		editor.Size(editor.Right(), editor.Bottom()+1);
	}
}

void Ted::CmToolbar()
{
	if (ToolbarCmd.isToggled())	{
		SetAttribute(TOOLBAR);
		AdjustBorders();
		editor.Size(editor.Right(), editor.Bottom()-2);
		editor.Move(editor.Left(), editor.Top()+2);
		toolbar.Enable();
		toolbar.Show();
		Paint();
	}
	else	{
		ClearAttribute(TOOLBAR);
		AdjustBorders();
		toolbar.Hide();
		toolbar.Disable();
		editor.Move(editor.Left(), editor.Top()-2);
		editor.Size(editor.Right(), editor.Bottom()+2);
		Paint();
	}
}

// ---- Options/Insert Menu Command
void Ted::CmInsert()
{
    editor.SetInsertMode(InsertCmd.isToggled());
}

// ---- Options/Tabs Menu Command
void Ted::CmTabs()
{
	char tbm[4];
	ostrstream ost(tbm, 4);
	ost << min(MaxTab, editor.Tabs()) << '\0';
	String tb(tbm);
	if (InputBox(2, "Tabs:", tb, 2))	{
		int tbs = atoi(tb);
		if (tbs <= MaxTab)	{
			editor.SetTabs(tbs);
			Paint();
		}
	}
}

// ----- reset the editor focus when the application moves
void Ted::Move(int x, int y)
{
    Application::Move(x, y);
    editor.SetFocus();
}

// ----- resize the editor when the application resizes
void Ted::Size(int x, int y)
{
    int wd = x - Left();
    int ht = y - Top();
	wd = max(wd, MinimumWidth);
	ht = max(ht, MinimumHeight);
	x = Left()+wd;
	y = Top()+ht;

    editor.Hide();
    editor.Size(editor.Right()+(x-Right()),
                editor.Bottom()+(y-Bottom()));
    Application::Size(x, y);
    editor.SetFocus();
}

// ---- test for changes to the document before discarding
void Ted::TestChanged()
{
    if (editor.Changed())    {
        String msg(fname + " has changed. Save?");
        if (YesNo(msg))
            CmSave();
    }
}

// ---- test for changes before closing
void Ted::CloseWindow()
{
    TestChanged();
    Application::CloseWindow();
}

void Ted::Paint()
{
	if (ruler)	{
		int wd = editor.Width();
		int tb = editor.Tabs();
		String rul(wd, '-');
		for (int i = tb+1; i < wd; i += tb)
			rul[i] = '|';
		rul[0] = rul[wd-1] = ' ';
		WriteClientString(rul, 0, 0, HighlightFG(), HighlightBG());
	}
}

void Ted::Execute()
{
    while (desktop.DispatchEvents())	{
		static int col = -1;
		static int row = -1;
		if (col != editor.GetColumn() || row != editor.GetRow())	{
			col = editor.GetColumn();
			row = editor.GetRow();
			Border();
		}
	}
}

void Ted::Border()
{
	Application::Border();
	char msg[30];
	ostrstream ost(msg, 30);
	ost << "Line: " << editor.GetRow() << " Col: " << editor.GetColumn() << '\0';
	StatusMessage(msg);
}

