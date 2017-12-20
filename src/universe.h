// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <vector>
#include <string>
#include <utility>
#include "mesh.h"
#include "field.h"
#include "jacobian.h"
#include "memory.h"
#include "operation.h"
#include "densematrix.h"
#include "selector.h"
#include "hierarchicalformfunction.h"
#include "hierarchicalformfunctioncontainer.h"

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
        
        // Store all !HIERARCHICAL! form function values that can be reused.
        // 'computedformfuncs[i].first' gives the ith form function type name.
        // 'computedformfuncs[i].second' gives a vector detailed below.
        // 'computedformfuncs[i].second[elemtypenum][0].first' gives the order
        // to which the form function has been interpolated and .second 
        // gives the interpolated values.
        static std::vector<std::pair< std::string, std::vector<std::vector< std::pair<int,hierarchicalformfunctioncontainer> >> >> computedformfuncs;
        // This function returns the requested form function value. In case
        // 'isreuseallowed' is true it reuses any already computed value.
        static hierarchicalformfunctioncontainer interpolateformfunction(std::string fftypename, int elementtypenumber, int interpolorder, std::vector<double> evaluationcoordinates);
        
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
        
        // If set to true the Gauss points weights product is not performed when assembling a formulation:
        static bool skipgausspointweightproduct;

        // If set to true the individual right handside contribution addresses and values are stored in 'rhsterms' when generating a formulation:
        static bool keeptrackofrhsassembly;
        // Every row in a given (int)densematrix corresponds to a given mesh element and every column to a shape function.
        // Do not forget to clear 'rhsterms' when you don't want to keep track anymore!
        static std::vector<std::pair<intdensematrix, densematrix>> rhsterms;

};

#endif
