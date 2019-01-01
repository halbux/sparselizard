// sparselizard - Copyright (C) 2017- A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef VECTORFIELDSELECT_H
#define VECTORFIELDSELECT_H

#include "memory.h"
#include "rawfield.h"
#include "rawvec.h"

class rawfield;
class rawvec;

class vectorfieldselect
{

	private:
       
        shared_ptr<rawvec> myrawvec;
        shared_ptr<rawfield> myrawfield;
        
	public:
        
        vectorfieldselect(shared_ptr<rawvec> vecin, shared_ptr<rawfield> fieldin) { myrawvec = vecin; myrawfield = fieldin; };
        
        shared_ptr<rawvec> getrawvector(void) { return myrawvec; };
        shared_ptr<rawfield> getrawfield(void) { return myrawfield; };
        
};

#endif
