// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object is a wrapper of the actual port object 'rawport' pointed
// by 'rawportptr' and to which most of the functions are redirected.

#ifndef PORT_H
#define PORT_H

#include <memory>
#include "rawport.h"

class rawport;

class port
{

    private:
    
        std::shared_ptr<rawport> rawportptr = NULL;
    
    public:
    
        port(void);
        // Provide the harmonics for a multiharmonic port:
        port(std::vector<int> harmonicnumbers);
        port(std::shared_ptr<rawport> rp);
        
        void setvalue(double portval);
        double getvalue(void);
        
        void setname(std::string name);
        std::string getname(void);
        
        std::vector<int> getharmonics(void);
        
        port harmonic(int harmonicnumber);
        port harmonic(std::vector<int> harmonicnumbers);
        port sin(int freqindex);
        port cos(int freqindex);
    
        std::shared_ptr<rawport> getpointer(void);
        
        void print(void);
    
        // Defining the +, -, * and / operators:
        expression operator+(void);
        expression operator-(void);

        expression operator+(port);
        expression operator-(port);
        expression operator*(port);
        expression operator/(port);

        expression operator+(double);
        expression operator-(double);
        expression operator*(double);
        expression operator/(double);

};

// Define the left version of the operators based on the right one.
expression operator+(double, port);
expression operator-(double, port);
expression operator*(double, port);
expression operator/(double, port);

#endif

