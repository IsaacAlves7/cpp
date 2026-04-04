// ------------- statbar.cpp

#include <time.h>
#include "statbar.h"

static Color col = {
	BLACK,CYAN,
	BLACK,CYAN,
	BLACK,CYAN,
	BLACK,CYAN
};

StatusBar::StatusBar(DFWindow *par)	: TextBox( par->ClientLeft(),
												par->Bottom(),
						  						1,
						  						par->ClientWidth(),
						  						par)
{
	windowtype = StatusbarWindow;
	SetAttribute(FRAMEWND);
	SetColor(col);
}

void StatusBar::ParentSized(int xdif, int ydif)
{
	Size(Right()+xdif, Top());
	Move(Left(), Bottom()+ydif);
}

void StatusBar::ClockTick()
{
    static Bool flipflop = False;
    time_t t = time(0);
    struct tm *now = localtime(&t);
    int hr = now->tm_hour > 12 ?
            now->tm_hour - 12 :
            now->tm_hour;
    if (hr == 0)
        hr = 12;
	int mn = now->tm_min;
	char hm[3];
	itoa(hr, hm, 10);
    String timestr(hm);
	timestr += ":";
	if (mn < 10)
		timestr += "0";
	itoa(mn, hm, 10);
	timestr += hm;
	timestr += (now->tm_hour > 11 ? "pm " : "am ");
	timestr += " ";
    /* ------- blink the : at one-second intervals ----- */
    if (flipflop)	{
        timestr[hr < 10 ? 1 : 2] = ' ';
	    flipflop = False;
	}
	else
	    flipflop = True;
	int x = ClientWidth() - 8;
	WriteClientString(timestr, x, 0, ClientFG(), ClientBG());
}

void StatusBar::StatusMessage(const String& Msg)
{
	int ln = Msg.Strlen();
	int x = (ClientWidth() - ln) / 2;
	WriteClientString(String(x, ' '), 0, 0, ClientFG(), ClientBG());
	WriteClientString(Msg, x, 0, ClientFG(), ClientBG());
	WriteClientString(String(x-8, ' '), x+ln, 0, ClientFG(), ClientBG());
}



