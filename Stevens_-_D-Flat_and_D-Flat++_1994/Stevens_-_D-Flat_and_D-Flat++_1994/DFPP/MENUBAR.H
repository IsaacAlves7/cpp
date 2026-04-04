// -------- menubar.h 

#ifndef MENUBAR_H
#define MENUBAR_H

#include "textbox.h"
#include "popdown.h"

class MenuBarItem    {
public:
    String *title;       // menu bar selection label
    int x1;              // 1st label position on bar
    int x2;              // last  "      "     "   "
    MenuSelection **ms;  // popdown selection list
    PopDown *popdown;    // popdown window
    void (DFWindow::*menuprep)();  // menu prep function
    MenuBarItem(const char *Title, MenuSelection **Ms = 0,
                            void (DFWindow::*MenuPrep)() = 0);
    ~MenuBarItem() { if (title) delete title; }
};

class MenuBar : public TextBox    {
    MenuBarItem *menuitems; // list of popdowns
    int menucount;          // count of popdowns
    int selection;          // current selection on the bar
    Bool ispoppeddown;      // True = a menu is down
    DFWindow *oldfocus;     // previous focus
    void Select();
    Bool AcceleratorKey(int key);
    Bool ShortcutKey(int key);
    void ParentSized(int xdif, int ydif);
public:
    MenuBar(MenuBarItem *MenuItems, DFWindow *par);
    ~MenuBar();
    // -------- menubar API messages
    void Keyboard(int key);
    void LeftButton(int mx, int my);
    Bool SetFocus();
    void ResetFocus();
    void Paint();
    void Select(int sel);
    void SetSelection(int sel);
};

#endif



