// -------- frame.h 

#ifndef FRAME_H
#define FRAME_H

#include "dfwindow.h"

class Frame : public DFWindow	{
	DFWindow *Host;		// window being moved or sized
	int diff;			// difference between mouse x and window corner
	int minx, maxx;		// x limits
	int miny, maxy;		// y limits
	Bool moving;		// True = moving, False = sizing
	void Paint() {}
	void LeftButton(int mx, int my);
	void ButtonReleased(int mx, int my);
	void MouseMoved(int mx, int my);
	void Show();
	void Hide();
	virtual void Keyboard(int key);
	void CloseFrame(Bool doOp);
public:
	Frame();
	void OpenFrame(DFWindow *Wnd, Bool Moving, int mx);
};

#endif




