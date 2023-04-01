// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef LOGS_H
#define LOGS_H

#include <sstream>
#include <iostream>
#include <string>

class logs
{

    private:

        std::ostringstream message;
    
    public:

        std::ostream& msg(void);
        
        void error(void);
    
};

#endif

