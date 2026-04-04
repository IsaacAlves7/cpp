// -------- dialog.h 

#ifndef DIALOG_H
#define DIALOG_H

#include "dfwindow.h"
#include "desktop.h"
#include "control.h"

class Control;

class Dialog : public DFWindow    {
    Bool isRunning;
    Bool okexit;
	DFWindow *ctlfocus;
    void OpenWindow();
    friend Control;
    void TestShortcut(int key);
    virtual Bool SetFocus();
protected:
    Dialog(const String& ttl, int lf, int tp, int ht, int wd,
            DFWindow *par = (DFWindow *)desktop.ApplWnd())
                : DFWindow(ttl, lf, tp, ht, wd, par)
            { OpenWindow(); }
    Dialog(const String& ttl, int ht, int wd,
            DFWindow *par = (DFWindow *)desktop.ApplWnd())
                : DFWindow(ttl, ht, wd, par)
            { OpenWindow(); }
    virtual void CloseWindow();
    virtual void Keyboard(int key);
    virtual void OKFunction();
    virtual void CancelFunction();
    virtual void HelpFunction();
    virtual void EnterFocus(DFWindow *Wnd);
public:
    virtual void Execute();
    Bool OKExit() { return okexit; }
};

#endif


