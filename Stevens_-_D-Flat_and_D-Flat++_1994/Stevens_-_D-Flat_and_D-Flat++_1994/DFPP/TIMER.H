// ----------- timer.h

#ifndef TIMER_H
#define TIMER_H

#include "dflatdef.h"

const int TIMER = 8;		// timer interrupt vector
const int MAXTIMERS = 5;	// maximum timers at one time

class Timer	{
	int timer;
	static Timer *timers[MAXTIMERS];
	static int timerct;
	static void interrupt (*OldTimer)(...);
	friend void interrupt NewTimer(...);
public:
	Timer();
	~Timer();
	Bool TimedOut()			 { return (Bool) (timer == 0); }
	void SetTimer(int ticks) { timer = ticks; }
	void DisableTimer() 	 { timer = -1; }
	Bool TimerRunning() 	 { return (Bool) (timer > 0); }
	void Countdown() 		 { --timer; }
	Bool TimerDisabled() 	 { return (Bool) (timer == -1); }
};

#endif


