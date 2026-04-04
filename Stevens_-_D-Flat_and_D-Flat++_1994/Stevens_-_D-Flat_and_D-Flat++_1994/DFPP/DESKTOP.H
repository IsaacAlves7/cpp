// ---------- desktop.h

#ifndef DESKTOP_H
#define DESKTOP_H

#include <dir.h>
#include "screen.h"
#include "cursor.h"
#include "keyboard.h"
#include "mouse.h"
#include "speaker.h"
#include "clock.h"
#include "clipbord.h"
#include "frame.h"

class DFWindow;
class TextBox;
class Application;
class Dialog;
class ErrorMessageBox;

const int ReserveSize = 512; // heap to reserve for exhaustion shutdown

const int CRIT = 0x24;
// ------- the interrupt function registers
struct Iregs {
	int bp,di,si,ds,es,dx,cx,bx,ax,ip,cs,fl;
};

class DeskTop {
	Application *apwnd;     // application window
	DFWindow *infocus;      // current window with the focus
	DFWindow *focuscapture; // current window with captured focus
	int drive;				// drive upon open
	char *path;				// subdirectory upon open
	unsigned errdrive;		// = drive if critical error occurs
	DFWindow *inWindow(DFWindow *wnd, int x, int y);
	// ------- the desktop devices
	Screen    sysscreen;    // the system screen
	Mouse     sysmouse;     // the system mouse
	Keyboard  syskeyboard;  // the system keyboard
	Cursor    syscursor;    // the system cursor
	Clock     sysclock;     // the system clock
	Speaker   sysspeaker;   // the system speaker
	Clipboard sysclipboard;	// the system clipboard
	Frame     sysframe;		// dummy for moving and sizing
	// ------ heap exhaustion support
	ErrorMessageBox *heaperror;		// used in new_handler
	char *reserve;			// held in reserve for new_handler
	void HeapError();		// called by new_handler
	Bool inError;			// semaphore for new_handler
	friend void NoMemory();
public:
	void SetHeapReserve(int rsv);

	DeskTop();
	~DeskTop();
	Application *ApplWnd()
		{ return apwnd; }
	void SetApplication(Application *ApWnd)
		{ apwnd = ApWnd; }
	Bool DispatchEvents();
	// --- focus methods
	DFWindow *InFocus()
		{ return infocus; }
	DFWindow *FocusCapture()
		{ return focuscapture; }
	void SetFocus(DFWindow *wnd)
		{ infocus = wnd; }
	void SetFocusCapture(DFWindow *wnd)
		{ focuscapture = wnd; }
	DFWindow *inWindow(int x, int y)
		{ return inWindow((DFWindow *)apwnd, x, y); }
	// --------- disk system methods
	void SetCriticalError(unsigned drive)
		{ errdrive = drive; }
	int TestCriticalError();
	String CurrentPath();
	int FindFirst(const String& path, ffblk *ff, int attrib);
	Bool SetDisk(int drive);
	// ------- the desktop devices
	Mouse&     mouse()      { return sysmouse;    }
	Screen&    screen()     { return sysscreen;   }
	Keyboard&  keyboard()   { return syskeyboard; }
	Cursor&    cursor()     { return syscursor;   }
	Clock&     clock()      { return sysclock;    }
	Speaker&   speaker()    { return sysspeaker;  }
	Frame&     frame()      { return sysframe; }
	Clipboard& clipboard()	{ return sysclipboard; }
};

extern DeskTop desktop;

#endif

