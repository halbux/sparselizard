// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef MYSTRING_H
#define MYSTRING_H

#include <iostream>
#include <vector>
#include <string>
#include <cctype>

class mystring
{

    private:

        int stringindex;
        std::string stringtoprocess;
    
    public:

        mystring(std::string);            
        // 'setstringindex' is optional. The index is 0 by default.
        void setstringindex(int stringindex);    

        std::string getstring(void);
        int getstringindex(void);            
        
        // 'getstringtonextwhitespace' gets the string between the current
        // string index and the next white space excluded. 'stringindex' 
        // is updated to the position of the next white space (or string end)
        // plus 1. If no white space is found then the string is returned 
        // till the end.
        std::string getstringtonextwhitespace(void);
        // 'jumptonextwhitespace' moves 'stringindex' to the position of the
        // next white space (or string end if none) plus 1. 
        void jumptonextwhitespace(void);
        
        // Similar to 'getstringtonextwhitespace'.
        std::string getstringtonextcomma(void);
    
        std::string getstringwhileletter(void);
    
};

#endif

