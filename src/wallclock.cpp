#include "wallclock.h"

wallclock::wallclock(void)
{
    starttime = std::chrono::steady_clock::now();
}

void wallclock::tic(void)
{
    accumulator = 0;
    ispaused = false;
    starttime = std::chrono::steady_clock::now();
}

double wallclock::toc(void)
{
    if (ispaused == false)
        return accumulator + std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - starttime).count();
    else
        return accumulator;
}

void wallclock::print(std::string toprint)
{
    if (toprint.size() > 0)
        std::cout << toprint << std::endl;

    double elapsedtime = toc();
    
    // If in the nanoseconds range:
    if (elapsedtime < 1000)
    {
        std::cout << elapsedtime << " ns" << std::endl;
        return;
    }
    elapsedtime = elapsedtime/1000;
    // If in the microseconds range:
    if (elapsedtime < 1000)
    {
        std::cout << elapsedtime << " us" << std::endl;
        return;
    }
    elapsedtime = elapsedtime/1000;
    // If in the milliseconds range:
    if (elapsedtime < 1000)
    {
        std::cout  << elapsedtime << " ms" << std::endl;
        return;
    }
    elapsedtime = elapsedtime/1000;
    // Otherwise we are in the seconds range:
    std::cout << elapsedtime << " s" << std::endl;
}

void wallclock::pause(void)
{
    if (ispaused == false)
    {
        accumulator += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - starttime).count();
        ispaused = true;
    }
}

void wallclock::resume(void)
{
    if (ispaused == true)
    {
        ispaused = false;
        starttime = std::chrono::steady_clock::now();
    }
}


