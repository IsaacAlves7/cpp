// --------- sound.h

class Sound	{
	void Wait(int n);
public:
	void Beep();
};

const int FREQUENCY = 100;
const long COUNT = 1193280L / FREQUENCY;
