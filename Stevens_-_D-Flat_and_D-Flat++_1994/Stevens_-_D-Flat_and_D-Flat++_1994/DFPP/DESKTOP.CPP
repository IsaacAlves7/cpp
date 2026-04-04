// --------------- desktop.cpp

#include <new.h>
#include <direct.h>
#include <ctype.h>
#include "desktop.h"
#include "applicat.h"
#include "notice.h"

DeskTop desktop;

// ------ critical error interrupt service routine
void interrupt newcrit(Iregs ir)
{
    if (!(ir.ax & 0x8000))
		desktop.SetCriticalError((ir.ax & 0xff) + 'A');
	ir.ax = 0;
}

void NoMemory()
{
	desktop.HeapError();
	exit(-1);
}

DeskTop::DeskTop()
{
	apwnd = 0;
	infocus = 0;
	focuscapture = 0;
	set_new_handler(NoMemory);
	syscursor.Hide();
	drive = getdisk();
	path = new char[129];
	getcwd(path, 129);
	errdrive = 0;
	setvect(CRIT, (void interrupt (*)(...))newcrit);
	heaperror = new ErrorMessageBox("Heap exhausted", "D-Flat++");
	reserve = new char[ReserveSize];
	inError = False;
}

DeskTop::~DeskTop()
{
	syscursor.Show();
	setdisk(drive);
	chdir(path);
	delete [] path;
	delete heaperror;
	delete reserve;
}

void DeskTop::SetHeapReserve(int rsv)
{
	delete reserve;
	reserve = new char[rsv];
}

Bool DeskTop::DispatchEvents()
{
	syskeyboard.DispatchEvent();
	sysmouse.DispatchEvent();
	sysclock.DispatchEvent();
	return (Bool) (apwnd != 0);
}

void DeskTop::HeapError()
{
	if (!inError)	{
		inError = True;
		delete [] reserve;
		reserve = 0;
		heaperror->Execute();
		if (apwnd != 0)
			apwnd->CloseWindow();
	}
}

// -------- test for critical errors
int DeskTop::TestCriticalError()
{
	if (errdrive)    {
		String ermsg(String("Drive ") + String(1, errdrive) + String(": Retry?"));
		errdrive = 0;
		if (YesNo(ermsg, "Disk Error"))
			return 1;
		else
			return 2;
	}
	return 0;
}

int DeskTop::FindFirst(const String& path, ffblk *ff, int attrib)
{
	int ax, err;
    do	{
        ax = findfirst(path, ff, attrib);
    	err = TestCriticalError();
	} while (err == 1);
	if (err)
		ax = -1;
    return ax;
}

String DeskTop::CurrentPath()
{
    char path[129];
	int err;
    do	{
	    _getdcwd(0, path, 129);
    	err = TestCriticalError();
	} while (err == 1);
	if (err)
		*path = '\0';
    return String(path);
}

Bool DeskTop::SetDisk(int drive)
{
    int currdrive = getdisk();
	int newdrive = toupper(drive) - 'A';
    setdisk(newdrive);
	// --- test to see if there is such a disk
	if (newdrive != getdisk())	{
		ErrorMessage("No such drive");
		return False;
	}
	// --- test to see if disk is accessible
	String cp = CurrentPath();
	if (cp.Strlen() == 0)	{
		setdisk(currdrive);
		return False;
	}
	return True;
}

inline Bool inShadow(Rect& rc, int x, int y)
{
	return (Bool)
	(((x == rc.Right()+1) && (y > rc.Top() && y < rc.Bottom()+2))
							||
	((y == rc.Bottom()+1) && (x > rc.Left() && x < rc.Right()+2)));
}

// ---- find the window that coordinates are in
DFWindow *DeskTop::inWindow(DFWindow *Wnd, int x, int y)
{
	DFWindow *hWnd;
	if (Wnd != 0)	{
		// ------ search window's children
		DFWindow *sWnd = Wnd->Last();
		while (sWnd != 0)	{
			hWnd = inWindow(sWnd, x, y);
			if (hWnd != 0)
				return hWnd;
			sWnd = sWnd->Prev();
		}
		// ---- test this window
		if (!Wnd->isVisible() || !Wnd->WindowRect().Inside(x, y))
			Wnd = 0;
	}
	return Wnd;
}


