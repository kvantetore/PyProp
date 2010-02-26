#ifndef TIMER_H
#define TIMER_H

class Timer
{
private:
	double time;

public:
	Timer() : time(0) {}

	void Start();
	void Stop();

	operator double() const
	{
		return time;
	}
};

#endif

