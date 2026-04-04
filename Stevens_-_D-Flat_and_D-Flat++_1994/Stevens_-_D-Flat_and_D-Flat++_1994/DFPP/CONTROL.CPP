// ------------ control.cpp

#include "control.h"
#include "desktop.h"

unsigned Control::LatestShortcut;

void Control::OpenControl()
{
	Enable();
	shortcut = LatestShortcut;
	LatestShortcut = 0;
	if (Parent())
		// --- x/y is relative to parent
		ChangePosition(Left() + Parent()->Left(),
						Top() + Parent()->Top());
}

void Control::Keyboard(int key)
{
    switch (key)    {
		case ESC:
			if (Parent())
				if (Parent()->WindowType() != DialogWindow)
					return;
			break;
        case '\t':
			key = ALT_F6;
			break;
        default:
            break;
    }
    DFWindow::Keyboard(key);
}

Bool Control::SetFocus()
{
	if (isEnabled())	{
		DFWindow::SetFocus();
		return True;
	}
	return False;
}

void Control::ShortcutSelect()
{
	SetFocus();
}

void Control::Select()
{
	if (Parent())
		Parent()->ControlSelected((DFWindow *)this);
}

void Control::Choose()
{
	if (Parent())
		Parent()->ControlChosen((DFWindow *)this);
}


