// ------------- focus.cpp

#include "desktop.h"
#include "dfwindow.h"

void DFWindow::CaptureFocus()
{
	if (this != desktop.FocusCapture())	{
		prevcapture = desktop.FocusCapture();
		prevfocus = desktop.InFocus();
		desktop.SetFocusCapture(this);
		if (this != desktop.InFocus())
			SetFocus();
	}
}

void DFWindow::ReleaseFocus()
{
	if (this == desktop.FocusCapture())	{
		desktop.SetFocusCapture(prevcapture);
		prevcapture = 0;
		if (prevfocus != 0)
			prevfocus->SetFocus();
		prevfocus = 0;
	}
}

Bool DFWindow::SetFocus()
{
	if (this != desktop.InFocus())	{
		// --- notify parent that child is taking the focus
		if (Parent() != 0)
			Parent()->EnterFocus(this);
		if (desktop.InFocus() != 0)
			desktop.InFocus()->ResetFocus();
		desktop.SetFocus(this);
		if (Parent() == 0 || Parent()->windowtype != DialogWindow)	{
			Dequeue();		// move the window to the
			Enqueue();		// top of the focus queue
		}
		Show();
	}
	return True;
}

void DFWindow::ResetFocus()
{
	desktop.SetFocus(0);
	Border();
	if (Parent() != 0)
		// --- notify parent that child lost the focus
		Parent()->LeaveFocus(this);
}

// --- set the focus to the next eligible sibling window
void DFWindow::NextSiblingFocus()
{
	DFWindow *Wnd = desktop.InFocus();
	while (Wnd != 0)	{
		if (Wnd->Next() == 0)	{
			if (Wnd->Parent() != 0)
				Wnd = Wnd->Parent()->First();
		}
		else
			Wnd = Wnd->Next();
		if (Wnd == desktop.InFocus())
			break;
		if (Wnd != 0)
			if (Wnd->SetFocus())
				break;
	}
}

// --- set the focus to the next previous eligible sibling window
void DFWindow::PrevSiblingFocus()
{
	DFWindow *Wnd = desktop.InFocus();
	while (Wnd != 0)	{
		if (Wnd->Prev() == 0)	{
			if (Wnd->Parent() != 0)
				Wnd = Wnd->Parent()->Last();
		}
		else
			Wnd = Wnd->Prev();
		if (Wnd == desktop.InFocus())
			break;
		if (Wnd != 0)
			if (Wnd->SetFocus())
				break;
	}
}

// --------- add the window to the linked list
void DFWindow::Enqueue()
{
	family.Append(this);
}

// --------- remove the window from the linked list
void DFWindow::Dequeue()
{
	family.Remove(this);
}




