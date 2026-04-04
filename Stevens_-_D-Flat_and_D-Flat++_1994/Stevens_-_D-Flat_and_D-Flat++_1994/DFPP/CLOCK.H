// --------- clock.h

#ifndef CLOCK_H
#define CLOCK_H

#include "timer.h"

class Clock    {
    Timer clocktimer;
public:
    Clock();
    void DispatchEvent();
};

#endif


