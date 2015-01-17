#ifndef L1_TIMERS
#define L1_TIMERS

#if defined _WIN32 || defined _WIN64

class Timer {
    public:
    Timer();
    unsigned long currentTicks();
    double currentSeconds();
    void start();
    void reset();

    __int64 timer_tickfreq;
    __int64 timer_start;
};

#else
        // support other platforms here
#include <sys/timeb.h>


class Timer {
    public:
    Timer();
    unsigned long currentTicks();
    double currentSeconds();
    void start();
    void reset();

    unsigned long timer_tickfreq;
    unsigned long timer_start;
};

#endif

#endif
