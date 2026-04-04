// ------- menusel.cpp

#include "menusel.h"

MenuSelection SelectionSeparator(SEPARATOR);

void MenuSelection::NullSelection()
{
    cmdfunction = 0;
    type = NORMAL;
    cascade = 0;
    isenabled = True;
    accelerator = 0;
    cascade = 0;
    toggle = Off;
}

void MenuSelection::CommonConstructor(
                               const String& Label,
                               int Accelerator,
			                   void (Application::*CmdFunction)(),
                               Bool Active,
                               MenuType Type,
                               Toggle Tgl,
                               MenuSelection **Cascaders)
{
    NullSelection();
    label = Label;
    accelerator = Accelerator;
    cmdfunction = CmdFunction;
    isenabled = Active;
    type = Type;
    toggle = Tgl;
    cascaders = Cascaders;
}

MenuSelection::MenuSelection(  const String& Label,
			                   void (Application::*CmdFunction)(),
                               int Accelerator,
                               Bool Active )
{
    CommonConstructor(Label, Accelerator, CmdFunction,
                                Active, NORMAL, Off);
}

MenuSelection::MenuSelection(  const String& Label,
			                   void (Application::*CmdFunction)(),
                               Toggle Tgl,
                               int Accelerator,
                               Bool Active)
{
    CommonConstructor(Label, Accelerator, CmdFunction,
                            Active, TOGGLE, Tgl);
}

MenuSelection::MenuSelection(const String& Label,
                            MenuSelection **Cascaders,
                            int Accelerator,
                            Bool Active )
{
    CommonConstructor(Label, Accelerator, 0,
                            Active, CASCADER, Off, Cascaders);
}

MenuSelection::MenuSelection(const String& Label,Toggle Tgl,
                             int Accelerator, Bool Active)
{
    CommonConstructor(Label, Accelerator, 0,
                            Active, TOGGLE, Tgl);
}

MenuSelection::MenuSelection(MenuType Type)
{
    NullSelection();
    type = Type;
}



