#include "timer.h"

void Timer::Start()
{
    gettimeofday(&start, nullptr);
}

long Timer::StopAndPrint(const std::string &label)
{
    struct timeval end;
    gettimeofday(&end, nullptr);

    long seconds = end.tv_sec - start.tv_sec;
    long useconds = end.tv_usec - start.tv_usec;
    long mtime = seconds * 1000000 + useconds;
    long mtime_sec = mtime / 1000000;
    float mtime_min = static_cast<float>(mtime_sec) / 60;

    printf("%s: %ld (%ld sec, %f min)\n", label.c_str(), mtime, mtime_sec, mtime_min);
    return mtime;
}

long Timer::StopTime()
{
    struct timeval end;
    gettimeofday(&end, nullptr);

    long seconds = end.tv_sec - start.tv_sec;
    long useconds = end.tv_usec - start.tv_usec;
    long mtime = seconds * 1000000 + useconds;
    return mtime;
}

void Timer::PrintTime(const std::string &label, long mtime)
{
    long mtime_sec = mtime / 1000000;
    float mtime_min = static_cast<float>(mtime_sec) / 60;
    printf("%s: %ld (%ld sec, %f min)\n", label.c_str(), mtime, mtime_sec, mtime_min);
}
