#include "timer.h"
#if defined _WIN32 || defined _WIN64
#include <windows.h>

Timer::Timer() {       
    timer_tickfreq = 1; // by default
    timer_start = 0;

    QueryPerformanceFrequency((_LARGE_INTEGER *)&timer_tickfreq);
}

void Timer::start() {
    QueryPerformanceCounter((_LARGE_INTEGER *)&timer_start);
}
void Timer::reset() { start(); }

unsigned long Timer::currentTicks() {
    __int64 elapsed = 0;
    QueryPerformanceCounter((_LARGE_INTEGER *)&elapsed);
    return (unsigned long)(elapsed - timer_start);
}

double Timer::currentSeconds() {
    return (timer_tickfreq != 0) ? (double)currentTicks()/(double)(timer_tickfreq) : 0;
}

#else

unsigned long GetMilliCount()
{
  // Something like GetTickCount but portable
  // It rolls over every ~ 12.1 days (0x100000/24/60/60)
  // Use GetMilliSpan to correct for rollover
  timeb tb;
  ftime( &tb );
  unsigned long nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}

unsigned long GetMilliSpan( unsigned long nTimeStart )
{
  int nSpan = (int) GetMilliCount() - nTimeStart;
  if ( nSpan < 0 )
    nSpan += 0x100000 * 1000;
  return nSpan;
}

Timer::Timer() {       
    timer_tickfreq = 1000; // by default
    timer_start = 0;

    //QueryPerformanceFrequency((_LARGE_INTEGER *)&timer_tickfreq);
}

void Timer::start() {
    //QueryPerformanceCounter((_LARGE_INTEGER *)&timer_start);
    timer_start = GetMilliCount();
}
void Timer::reset() { start(); }

unsigned long Timer::currentTicks() {
    //__int64 elapsed = 0;
    //QueryPerformanceCounter((_LARGE_INTEGER *)&elapsed);
    //return (unsigned long)(elapsed - timer_start);

    return GetMilliSpan(timer_start);
}

double Timer::currentSeconds() {
    return (timer_tickfreq != 0) ? (double)currentTicks()/(double)(timer_tickfreq) : 0;
}
#endif

