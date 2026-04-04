// -------- editbox.h 

#ifndef EDITBOX_H
#define EDITBOX_H

#include "textbox.h"

class EditBox : public TextBox	{
	void OpenWindow();
protected:
    int column;       // Current column
    Bool changed;     // True if text has changed
	Bool insertmode;
	virtual void Home();
	virtual void End();
	virtual void NextWord();
	virtual void PrevWord();
	virtual void Forward();
	virtual void Backward();
    virtual Bool PageRight();
    virtual Bool PageLeft();
	virtual void DeleteCharacter();
	virtual void InsertCharacter(int key);
	virtual void ExtendBlock(int x, int y);
public:
	EditBox(const String& ttl, int lf, int tp, int ht, int wd, DFWindow *par=0)
						: TextBox(ttl, lf, tp, ht, wd, par)
			{ OpenWindow(); }
	EditBox(const String& ttl, int ht, int wd, DFWindow *par=0)
						: TextBox(ttl, ht, wd, par)
			{ OpenWindow(); }
	EditBox(int lf, int tp, int ht, int wd, DFWindow *par=0)
						: TextBox(lf, tp, ht, wd, par)
			{ OpenWindow(); }
	EditBox(int ht, int wd, DFWindow *par=0) : TextBox(ht, wd, par)
			{ OpenWindow(); }
	EditBox(const String& ttl)	: TextBox(ttl)
			{ OpenWindow(); }
	// -------- API messages
	virtual Bool SetFocus();
	virtual void ResetFocus();
	virtual void SetCursor(int x, int y);
	virtual void SetLeftMargin(int x)
			{ wleft = x; }
	virtual Bool ResetCursor();
	virtual void SetCursorSize();
	virtual unsigned char CurrentChar()
			{ return (*text)[column]; }
	virtual unsigned CurrentCharPosition()
			{ return column; }
	virtual Bool AtBufferStart()
			{ return (Bool) (column == 0); }
	virtual void Keyboard(int key);
	virtual void Move(int x, int y);
	virtual void PaintCurrentLine()
			{ Paint(); }
	virtual void ClearText();
    virtual void LeftButton(int mx, int my);
	Bool Changed()
			{ return changed; }
	void ClearChanged()
			{ changed = False; }
	Bool InsertMode()
			{ return insertmode; }
	void SetInsertMode(Bool imode)
			{ insertmode = imode; ResetCursor(); }
	virtual int GetColumn() const { return column; }
	virtual void DeleteSelectedText();
	virtual int CurrentPosition() const
		{ return (int) (TextLine(GetRow()) - *text) + column; }
	virtual void InsertText(const String& txt);
	void Cut();
	void Copy();
	void Paste();
};

inline Bool isWhite(int ch)
{
	return (Bool) (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}

#endif


