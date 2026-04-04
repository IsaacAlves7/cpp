// -------- popdown.h 

#ifndef POPDOWN_H
#define POPDOWN_H

#include "desktop.h"
#include "listbox.h"

const unsigned char LEDGE          = '\xc3';
const unsigned char REDGE          = '\xb4';
const unsigned char CASCADEPOINTER = '\x10';

inline unsigned char CheckMark()
{
    return desktop.screen().Height() == 25 ? 251 : 4;
}

class MenuSelection;
class MenuBar;

class PopDown : public ListBox    {
    MenuSelection **selections; // array of selections
    Bool isopen;                // True = menu is open
    Bool iscascaded;            // True = menu is cascaded
    int menuwidth;              // width of menu
    int menuheight;             // height of menu

    void BuildMenuLine(int sel);
    void MenuDimensions();
    void DisplayMenuLine(int lno);
    Bool ShortcutKey(int key);
protected:
    void ClearSelection();
public:
    PopDown(DFWindow *par, MenuSelection **Selections = 0);
    // -------- listbox API messages
    void CloseWindow();
    void OpenMenu(int left, int top);
    void CloseMenu(Bool SendESC = False);
    void Show();
    void Paint();
    void Border();
    void Keyboard(int key);
    void ShiftChanged(int sk);
    void ButtonReleased(int mx, int my);
    void LeftButton(int mx, int my);
    void DoubleClick(int mx, int my);
    void Choose();
    void SetSelection(int sel);
    Bool isOpen() { return isopen; }
    Bool AcceleratorKey(int key);
    Bool ParentisMenu(DFWindow& wnd);
    Bool ParentisMenu() { return ParentisMenu(*this); }
};

#endif



