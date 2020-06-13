// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef WALLCLOCK_H
#define WALLCLOCK_H

#include <iostream>
#include <chrono>
#include <sstream>

class wallclock
{
    private:

        std::chrono::steady_clock::time_point starttime;
        
        double accumulator = 0;
        
        bool ispaused = false;
        
    public:
    
        // The constructor acts like 'tic':
        wallclock(void);
        // Reset clock:
        void tic(void);
        // Returns time since previous tic or object creation in nanoseconds.
        // Takes into account any accumulated time via pause/resume.
        double toc(void);
        // Prints the time since previous tic or object creation in an adapted unit.
        // Also prints the message in the argument string (if any).
        void print(std::string toprint = "");
        
        // Pause the clock:
        void pause(void);
        // Resume the clock:
        void resume(void);
};

#endif
