// ---------- control.h

#ifndef CONTROL_H
#define CONTROL_H

#include "dfwindow.h"

class TextBox;
class Dialog;

class Control : public DFWindow    {
    Bool enabled;       // true if control is enabled
protected:
	unsigned shortcut;
	static unsigned LatestShortcut;
public:
    Control(const String& ttl,int lf,int tp,int ht,int wd,DFWindow *par)
                        : DFWindow(ttl, lf, tp, ht, wd, par)
        { OpenControl(); }
    Control(const String& ttl, int ht, int wd, DFWindow *par)
                        : DFWindow(ttl, ht, wd, par)
        { OpenControl(); }
    Control(int lf, int tp, int ht, int wd, DFWindow *par)
                        : DFWindow(lf, tp, ht, wd, par)
        { OpenControl(); }
    Control(int ht,int wd,DFWindow *par) : DFWindow(ht,wd,par)
        { OpenControl(); }
    Control(const String& ttl) : DFWindow(ttl)
        { OpenControl(); }
    virtual void Keyboard(int key);
    void OpenControl();
    void Enable()      { enabled = True; }
    void Disable()     { enabled = False; }
    Bool isEnabled()   { return enabled; }
	virtual void ShortcutSelect();
    virtual Bool SetFocus();
	unsigned Shortcut() { return shortcut; }
	virtual void Select();
	virtual void Choose();
};

#endif

