// ------------- textbox.h

#ifndef TEXTBOX_H
#define TEXTBOX_H

#include "control.h"

const int SHORTCUTCHAR = '~';

class ScrollBar;

enum Marking { NotMarking, MouseMarking, KeyboardMarking };

class TextBox : public Control    {
    ScrollBar *hscrollbar;  // horizontal scroll bar
    ScrollBar *vscrollbar;  // vertical scroll bar
    unsigned *textpointers; // -> list of line offsets
    Bool resetscrollbox;
    void OpenWindow();
	Bool buttondown;
	void PadString(int lno, int datlen);
    virtual void ShiftChanged(int sk);
protected:
	Marking textmarking;
	int ancx, ancy;
    // ---- text buffer
    String *text;           // window text
    int wlines;             // number of lines of text
    int textwidth;          // width of longest line in textbox
	unsigned textlength;	// maximum allowed text length (0 = no max)
    // ---- text display
    int wtop;               // text line on top of display
    int wleft;              // left position in window viewport
    int BlkBegLine;         // beginning line of marked block
    int BlkBegCol;          // beginning column of marked block
    int BlkEndLine;         // ending line of marked block
    int BlkEndCol;          // ending column of marked block
    int shortcutfg;         // shortcut key color
	// ---- Selected Text
	void SetAnchor(int x, int y);
	virtual void ExtendBlock(int x, int y);
    const char *TextLine(int line) const
        { return (const char *)(*text) + *(textpointers+line); }
    int DisplayShortcutField(
        const String& sc, int x, int y, int fg, int bg);
    void WriteShortcutLine(int lno, int fg, int bg);
    void WriteTextLine(int lno);
    void WriteTextLine(int lno, int fg, int bg);
    void PaintCurrentLine()
        { WriteTextLine(GetRow()); }
    virtual void BuildTextPointers();
	unsigned *TextPointers() { return textpointers; }
    void SetScrollBoxes();
    virtual void MouseMoved(int mx, int my);
    virtual void ButtonReleased(int mx, int my);
	void KeepInText(int& mx, int& my);
	int LineLength(int lno);
	void TestMarking();
	virtual int GetRow() const { return 0; }
	virtual int GetColumn() const { return 0; }
public:
    TextBox(const String& ttl,int lf,int tp,int ht,int wd,DFWindow *par=0)
                        : Control(ttl, lf, tp, ht, wd, par)
            { OpenWindow(); }
    TextBox(const String& ttl, int ht, int wd, DFWindow *par=0)
                        : Control(ttl, ht, wd, par)
            { OpenWindow(); }
    TextBox(int lf, int tp, int ht, int wd, DFWindow *par=0)
                        : Control(lf, tp, ht, wd, par)
            { OpenWindow(); }
    TextBox(int ht, int wd, DFWindow *par=0) : Control(ht,wd,par)
            { OpenWindow(); }
    TextBox(const String& ttl) : Control(ttl)
            { OpenWindow(); }
    // -------- textbox API messages
    virtual Bool ScrollUp();
    virtual Bool ScrollDown();
    virtual Bool ScrollRight();
    virtual Bool ScrollLeft();
    virtual Bool PageUp();
    virtual Bool PageDown();
    virtual Bool PageRight();
    virtual Bool PageLeft();
    virtual void Home();
    virtual void End();
    virtual void CloseWindow();
    virtual void AddText(const String& txt);
    virtual void SetText(const String& txt);
    virtual void SetTextLength(unsigned int len);
	virtual const String GetText() { return *text; }
    virtual void ClearText();
    virtual void Show();
    virtual void Paint();
    virtual void Keyboard(int key);
    virtual void LeftButton(int mx, int my);

    void ExtractTextLine(String& ln, int lno);
	unsigned int TextLength() const
		{ return text ? text->Strlen() : 0; }
	void SelectText(int lf, int tp, int rt, int bt);
	virtual String GetSelectedText() const;
	virtual void DeleteSelectedText();
	int BlockBeginOffset() const;
	int BlockEndOffset() const;
	int BlockLength() const
		{ return BlockEndOffset() - BlockBeginOffset(); }
    void ClearTextBlock();
	Bool TextBlockMarked() const
		{ return (Bool)
			!(BlkBegLine == BlkEndLine && BlkBegCol == BlkEndCol); }
    void HorizontalPagePosition(int pct);
    void VerticalPagePosition(int pct);
	int LineCount() const { return wlines; }
	int TextWidth() const { return textwidth; }
	virtual void SetScrollBars();
};

#endif

