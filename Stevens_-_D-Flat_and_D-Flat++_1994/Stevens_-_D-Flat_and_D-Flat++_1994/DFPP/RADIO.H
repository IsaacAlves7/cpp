// ----------- radio.h

#ifndef RADIO_H
#define RADIO_H

#include "button.h"

class RadioButton : public Button    {
public:
    RadioButton(const String& lbl, int lf, int tp, DFWindow *par);
    // -------- API messages
    virtual void InvertButton();
    virtual void PushButton();
};

#endif

