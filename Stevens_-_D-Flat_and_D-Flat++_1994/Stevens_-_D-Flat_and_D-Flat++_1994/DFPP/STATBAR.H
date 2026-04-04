// -------- statbar.h 

#ifndef STATBAR_H
#define STATBAR_H

#include "textbox.h"

class StatusBar : public TextBox	{
	void ParentSized(int xdif, int ydif);
protected:		 
public:
	StatusBar(DFWindow *par);
	void ClockTick();
	void StatusMessage(const String& Msg);
	Bool SetFocus() { return False; }
};

#endif



