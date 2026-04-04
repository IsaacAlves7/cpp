// -------- timer.cpp

#include <dos.h>
#include "timer.h"

void interrupt (*Timer::OldTimer)(...);

Timer *Timer::timers[MAXTIMERS];
int Timer::timerct = 0;

// ------- timer interrupt service routine
void interrupt NewTimer(...)
{
	// --- countdown all running timers
	for (int i = 0; i < MAXTIMERS; i++)
		if (Timer::timers[i] != 0)
		    if (Timer::timers[i]->TimerRunning())
    		    Timer::timers[i]->Countdown();
	// --- chain to the old interrupt vector
    (*Timer::OldTimer)();
}

Timer::Timer()
{
	timer = -1;
	// --- if first timer, hook and chain the interrupt vector
	if (timerct == 0)	{
	    OldTimer = getvect(TIMER);
    	setvect(TIMER, NewTimer);
	}
	// --- add the timer to the table
	timers[timerct++] = this;
}

Timer::~Timer()
{
	// --- remove the timer from the table
	for (int i = 0; i < MAXTIMERS; i++)	{
		if (this == timers[i])	{
			timers[i] = 0;
			break;
		}
	}
	// --- if last one, restore the interrupt vector
	if (--timerct == 0)
		setvect(TIMER, OldTimer);
}


