#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <sys/time.h>

class timer{
public:
  timer(void);
  void start(void);
  void stop(void);
  double time(void);
  
private:
  timeval start_time;
  double run_time;
  bool running;
};

#endif