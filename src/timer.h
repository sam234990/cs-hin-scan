#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>
#include <cstdio>
#include <string>

class Timer {
public:
    void Start();
    long StopAndPrint(const std::string& label); 
    long StopTime();
    static void PrintTime(const std::string& label, long mtime);
private:
    struct timeval start;
};

#endif // TIMER_H
