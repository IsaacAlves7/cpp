// ------- clock.cpp

#include "desktop.h"
#include "applicat.h"

Clock::Clock()
{
    clocktimer.SetTimer(0);
}

void Clock::DispatchEvent()
{
    if (clocktimer.TimedOut())    {
        // -------- reset the timer
        clocktimer.SetTimer(19);    // approx. 19 tics/second
        // -------- post the clock event
        if (desktop.ApplWnd() != 0)
            desktop.ApplWnd()->ClockTick();
    }
}

