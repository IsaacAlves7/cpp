// ------------- radio.cpp

#include "radio.h"
#include "desktop.h"

RadioButton::RadioButton(const String& lbl, int lf, int tp,
        DFWindow *par) : Button(lbl, lf, tp, par)
{
    windowtype = RadioButtonWindow;
    setchar = 7;
}

void RadioButton::InvertButton()
{
    PushButton();
}

void RadioButton::PushButton()
{
    int ht = desktop.screen().Height();
    DFWindow **rd = new DFWindow *[ht];
    for (int i = 0; i < ht; i++)
        rd[i] = 0;
    // - build a table of radio buttons at the same x coordinate
    DFWindow *sib = Parent()->First();
    while (sib != 0)    {
        if (sib->WindowType() == RadioButtonWindow)    {
            if (sib->Left() == Left())    {
                int tp = sib->Top();
                if (tp < ht)
                    rd[tp] = sib;
            }
        }
        sib = sib->Next();
    }
    // ----- find the start of the radiobutton group
    i = Top();
    while (i >= 0 && rd[i] != 0)
        --i;
    // ---- ignore everthing before the group
    while (i >= 0)
        rd[i--] = 0;
    // ----- find the end of the radiobutton group
    i = Top();
    while (i < ht && rd[i] != 0)
        i++;
    // ---- ignore everthing past the group
    while (i < ht)
        rd[i++] = 0;
    // ------ release all the radio buttons in the group
    for (i = 0; i < ht; i++)
        if (rd[i] != 0)
            ((RadioButton *)rd[i])->ReleaseButton();
    delete [] rd;
    // ----- set the chosen radio button
    Button::PushButton();
}
