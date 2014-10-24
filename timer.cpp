#include "timer.h"

static timeval getTime(void){
  timeval tv;
  gettimeofday(&tv, NULL);
  return tv;
}

static double timeDiff(timeval end, timeval start) {
  return (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
}

timer::timer(void) : run_time(0), running(false){
}

void timer::start(void){
  if (running == true) return;
  start_time = getTime();
  running = true;
}

void timer::stop(void){
  if (running == false) return;
  timeval stop_time = getTime();
  run_time += timeDiff(stop_time, start_time);
  running = false;
}

double timer::time(void){
  if (running == true){
    timeval stop_time = getTime();
    return run_time + timeDiff(stop_time, start_time);
  } else
    return run_time;
}