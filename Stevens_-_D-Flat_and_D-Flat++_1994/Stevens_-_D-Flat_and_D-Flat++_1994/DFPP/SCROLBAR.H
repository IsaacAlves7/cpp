// -------- scrolbar.h 

#ifndef SCROLBAR_H
#define SCROLBAR_H

#include "dfwindow.h"

class TextBox;

/* ------------- scroll bar characters ------------ */
const unsigned char UPSCROLLBOX    = '\x1e';
const unsigned char DOWNSCROLLBOX  = '\x1f';
const unsigned char LEFTSCROLLBOX  = '\x11';
const unsigned char RIGHTSCROLLBOX = '\x10';
const unsigned char SCROLLBARCHAR  = 176 ;
const unsigned char SCROLLBOXCHAR  = 178;

enum BarPlane { HORIZONTAL, VERTICAL };

class ScrollBar : public DFWindow	{
	int scrollbox;			// position of scroll box
	Bool sliding;			// True while siding the scroll box
	BarPlane plane;			// HORIZONTAL or VERTICAL
	void ParentSized(int xdif, int ydif);
protected:
public:
	ScrollBar(BarPlane Plane, TextBox *par);
	// -------- ScrollBar API messages
	void Paint();
	Bool SetFocus();
	void LeftButton(int mx, int my);
	void MouseMoved(int mx, int my);
	void ButtonReleased(int mx, int my);
	void TextPosition(int TxPct);
	void MoveScrollBox(int sb);
};

#endif



