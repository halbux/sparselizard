// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef VECTORFIELDSELECT_H
#define VECTORFIELDSELECT_H

#include "memory.h"
#include "rawfield.h"
#include "rawvec.h"
#include "field.h"

class rawfield;
class rawvec;
class field;

class vectorfieldselect
{

    private:
       
        std::shared_ptr<rawvec> myrawvec;
        std::shared_ptr<rawfield> myrawfield;
        
    public:
        
        vectorfieldselect(std::shared_ptr<rawvec> vecin, std::shared_ptr<rawfield> fieldin) { myrawvec = vecin; myrawfield = fieldin; };
        
        std::shared_ptr<rawvec> getrawvector(void) { return myrawvec; };
        std::shared_ptr<rawfield> getrawfield(void) { return myrawfield; };
        
        void setdata(int physreg, field myfield, std::string op = "set");
        
};

#endif
