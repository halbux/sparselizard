// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <vector>
#include "mesh.h"
#include "field.h"
#include "jacobian.h"
#include "memory.h"
#include "operation.h"
#include "densematrix.h"

class mesh;
class jacobian;

class universe 
{
    private:
    
    public:

        static mesh* mymesh;
        
        static double currenttimestep;
        
        static double fundamentalfrequency;
        static double getfundamentalfrequency(void);
        
        // To allow reusing computed things:
        static bool isreuseallowed;
        static void allowreuse(void);
        // CLEANS::
        static void forbidreuse(void);
        
        static shared_ptr<jacobian> computedjacobian;
        
        // Store all operations that must be reused:
        static std::vector<shared_ptr<operation>> oppointers;
        static std::vector<shared_ptr<operation>> oppointersfft;
        // Store all computed values:
        static std::vector< std::vector<std::vector<densematrix>> > opcomputed;
        static std::vector< densematrix > opcomputedfft;
        
        // Returns -1 if not yet precomputed.
        static int getindexofprecomputedvalue(shared_ptr<operation> op);
        static int getindexofprecomputedvaluefft(shared_ptr<operation> op);
        // Returns a copy to avoid any modification of the data stored here:
        static std::vector<std::vector<densematrix>> getprecomputed(int index);
        static densematrix getprecomputedfft(int index);
        // Sets a copy to avoid any modification of the data stored here:
        static void setprecomputed(shared_ptr<operation> op, std::vector<std::vector<densematrix>> val);
        static void setprecomputedfft(shared_ptr<operation> op, densematrix val);

};

#endif
