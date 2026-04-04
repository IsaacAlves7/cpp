// ----------- button.h

#ifndef BUTTON_H
#define BUTTON_H

#include "textbox.h"

class Button : public TextBox    {
    Bool setting;
protected:
    char setchar;
    Button(const String& lbl, int lf, int tp, DFWindow *par);
public:
    // -------- API messages
    virtual Bool SetFocus();
    virtual void ResetFocus();
    virtual void Paint();
    virtual void Keyboard(int key);
    virtual void LeftButton(int mx, int my);
    virtual void InvertButton();
    virtual void PushButton();
    virtual void ReleaseButton();
    Bool Setting() { return setting; }
	virtual void ShortcutSelect() { InvertButton(); }
};

#endif

