// --------- speaker.h

#ifndef SPEAKER_H
#define SPEAKER_H

class Speaker    {
public:
    void Beep();
};

const int FREQUENCY = 100;
const long COUNT = 1193280L / FREQUENCY;

#endif


