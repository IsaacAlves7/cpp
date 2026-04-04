// ------- menusel.h

#ifndef MENUSEL_H
#define MENUSEL_H

#include <stdio.h>
#include "strings.h"
#include "dfwindow.h"

enum MenuType {
    NORMAL,
    TOGGLE,
    CASCADER,
    SEPARATOR
};

enum Toggle { Off, On };

class DFWindow;
class PopDown;

class MenuSelection    {
    String label;           // selection label
    void (Application::*cmdfunction)(); // selection function
    MenuType type;          // NORMAL, TOGGLE,
                            // CASCADER, SEPARATOR
    Bool isenabled;         // True = enabled, False = disabled
    int accelerator;        // accelerator key
    PopDown *cascade;       // cascaded menu window
    MenuSelection **cascaders; // cascaded menu selection list
    Toggle toggle;          // On or Off
    void NullSelection();   // Build a null selection
    friend PopDown;
    void CommonConstructor( const String& Label,
                            int Accelerator,
		                    void (Application::*CmdFunction)(),
                            Bool Active,
                            MenuType Type,
                            Toggle Tgl,
                            MenuSelection **Cascaders = 0);
public:
    MenuSelection(  const String& Label, 
                    void (Application::*CmdFunction)() = 0,
                    int Accelerator=0,
                    Bool Active=True );

    MenuSelection(  const String& Label,
                    void (Application::*CmdFunction)(),
                    Toggle Tgl,
                    int Accelerator=0,
                    Bool Active=True );

    MenuSelection(  const String& Label,
                    MenuSelection **Cascaders,
                    int Accelerator=0,
                    Bool Active=True );

    MenuSelection(  const String& Label,
                    Toggle Tgl,
                    int Accelerator=0,
                    Bool Active=True );

    MenuSelection(MenuType Type);

    Bool isEnabled()    { return isenabled; }
    void Enable()       { isenabled = True; }
    void Disable()      { isenabled = False; }
    void SetToggle()    { toggle = On; }
    void ClearToggle()  { toggle = Off; }
    void InvertToggle() { toggle = toggle == On ? Off : On; }
    Bool isToggled()    { return (Bool) (toggle == On); };
    Type()              { return type; }
};

extern MenuSelection SelectionSeparator;
extern MenuSelection SelectionTerminator;

#endif


