// -------- editor.h 

#ifndef EDITOR_H
#define EDITOR_H

#include "editbox.h"

const unsigned char Ptab = '\t'+0x80; // pseudo tab expansion
#define MaxTab 12                     // maximum tab width

class Editor : public EditBox    {
    int tabs;           // tab expansion value
    Bool wordwrapmode;  // True = wrap words
    int row;            // Current row
    void OpenWindow();
    void ScrollCursor();
    void AdjustCursorTabs(int key = 0);
    void ExtendBlock(int x, int y);
    void InsertTab();
    void AdjustTabInsert();
    void AdjustTabDelete();
    void WordWrap();
    Bool AtBufferStart()
        { return (Bool) (column == 0 && row == 0); }
protected:
    virtual void Upward();
    virtual void Downward();
    virtual void Forward();
    virtual void Backward();
    virtual void DeleteCharacter();
    virtual void BeginDocument();
    virtual void EndDocument();
    virtual Bool PageUp();
    virtual Bool PageDown();
    virtual void InsertCharacter(int key);
    virtual Bool ResetCursor();
    virtual void WriteString(const String& ln,
                        int x, int y, int fg, int bg);
    virtual void LeftButton(int mx, int my);
public:
    Editor(const String& ttl, int lf, int tp, int ht, int wd,
                                            DFWindow *par=0)
                        : EditBox(ttl, lf, tp, ht, wd, par)
            { OpenWindow(); }
    Editor(const String& ttl, int ht, int wd, DFWindow *par=0)
                        : EditBox(ttl, ht, wd, par)
            { OpenWindow(); }
    Editor(int lf, int tp, int ht, int wd, DFWindow *par=0)
                        : EditBox(lf, tp, ht, wd, par)
            { OpenWindow(); }
    Editor(int ht, int wd, DFWindow *par=0)
                        : EditBox(ht, wd, par)
            { OpenWindow(); }
    Editor(const String& ttl) : EditBox(ttl)
            { OpenWindow(); }
    virtual void Keyboard(int key);
    virtual unsigned char CurrentChar()
        { return *(TextLine(row) + column); }
    virtual unsigned CurrentCharPosition()
        { return (unsigned)
            ((const char *)
            (TextLine(row)+column) - (const char *) *text);
        }
    virtual void FormParagraph();
    virtual void AddText(const String& txt);
    virtual const String GetText();
    virtual void ClearText();
    virtual int GetRow() const { return row; }
    int Tabs()
        { return tabs; }
    void SetTabs(int t);
    Bool WordWrapMode()
        { return wordwrapmode; }
    void SetWordWrapMode(Bool wmode)
        { wordwrapmode = wmode; }
    virtual void DeleteSelectedText();
    virtual void InsertText(const String& txt);
	virtual void PageToLine(int y)
		{ wtop = y; }
	virtual void SetCursor(int x, int y);
};

#endif
