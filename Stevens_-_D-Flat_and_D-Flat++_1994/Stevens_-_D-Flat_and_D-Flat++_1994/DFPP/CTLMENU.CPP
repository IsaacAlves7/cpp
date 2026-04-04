// --------------- ctlmenu.cpp

#include "applicat.h"

#include "menusel.h"
#include "frame.h"

#define Ap void (Application::*)()

MenuSelection MoveCmd     ("~Move",     (Ap) &DFWindow::CtlMenuMove);
MenuSelection SizeCmd     ("~Size",     (Ap) &DFWindow::CtlMenuSize);
MenuSelection CloseDocCmd ("~Close  [Ctrl+F4]",
                                        (Ap) &DFWindow::CloseWindow, CTRL_F4);
MenuSelection CloseApCmd  ("~Close  [Alt+F4]",
                                        (Ap) &DFWindow::CloseWindow, ALT_F4);

MenuSelection *ControlMenu[8];

MenuBarItem CtlMenu[] = {
    MenuBarItem( "", ControlMenu ),
    MenuBarItem( 0, 0 )
};

void DFWindow::OpenCtlMenu()
{
    if (ctlmenu != 0)
        delete ctlmenu;

    int mn = 0;

    if (Attribute() & MOVEABLE)
        ControlMenu[mn++] = &MoveCmd;
    if (Attribute() & SIZEABLE)
        ControlMenu[mn++] = &SizeCmd;
    if (mn != 0)
        ControlMenu[mn++] = &SelectionSeparator;
    if (Parent() && WindowType() != DialogWindow)
        ControlMenu[mn++] = &CloseDocCmd;
    else
        ControlMenu[mn++] = &CloseApCmd;
    ControlMenu[mn] = 0;

    SizeCmd.Disable();
    MoveCmd.Disable();

    if (windowstate == OPEN)	{
        MoveCmd.Enable();
        SizeCmd.Enable();
    }
    ctlmenu = new PopDown(this, ControlMenu);
    ctlmenu->OpenMenu(Left()+1, Top()+1);
}

void DFWindow::DeleteCtlMenu()
{
    if (ctlmenu != 0)
        delete ctlmenu;
    ctlmenu = 0;
}

void DFWindow::CtlMenuMove()
{
    desktop.mouse().SetPosition(Left(), Top());
    desktop.frame().OpenFrame(this, True, Left());
}

void DFWindow::CtlMenuSize()
{
    desktop.mouse().SetPosition(Right(), Bottom());
    desktop.frame().OpenFrame(this, False, 0);
}


