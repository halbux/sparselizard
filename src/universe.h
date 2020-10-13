// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <vector>
#include <string>
#include <utility>
#include "rawmesh.h"
#include "field.h"
#include "jacobian.h"
#include "memory.h"
#include "operation.h"
#include "densematrix.h"
#include "selector.h"
#include "hierarchicalformfunction.h"
#include "hierarchicalformfunctioncontainer.h"
#include "vec.h"

class mesh;
class jacobian;

class universe 
{
    private:
    
    public:

        // 'Round-off noise level on the node coordinates:
        static double roundoffnoiselevel;
        
        static std::shared_ptr<rawmesh> mymesh;

        static bool isaxisymmetric;
        
        static double currenttimestep;
        
        static double fundamentalfrequency;
        static double getfundamentalfrequency(void);
        
        // Shift the physical region numbers by (physregdim+1) x physregshift when loading a mesh:
        static int physregshift;
        
        // Error estimators need evaluations across elements and can therefore not be efficiently evaluated in the usual way.
        // Error estimators will be updated if 'numallowedtimes' > 0 and their state number is different than the one here.
        static void allowestimatorupdate(bool allowitonce);
        static bool isestimatorupdateallowed(long long int statenumber);
        static long long int estimatorcalcstate;
        static int numallowedtimes;
        
        // To allow reusing computed things:
        static bool isreuseallowed;
        static void allowreuse(void);
        // CLEANS::
        static void forbidreuse(void);
        
        static std::shared_ptr<jacobian> computedjacobian;
        
        // Store all operations that must be reused:
        static std::vector<std::shared_ptr<operation>> oppointers;
        static std::vector<std::shared_ptr<operation>> oppointersfft;
        // Store all computed values:
        static std::vector< std::vector<std::vector<densematrix>> > opcomputed;
        static std::vector< densematrix > opcomputedfft;
        
        // Returns -1 if not yet precomputed.
        static int getindexofprecomputedvalue(std::shared_ptr<operation> op);
        static int getindexofprecomputedvaluefft(std::shared_ptr<operation> op);
        // Returns a copy to avoid any modification of the data stored here:
        static std::vector<std::vector<densematrix>> getprecomputed(int index);
        static densematrix getprecomputedfft(int index);
        // Sets a copy to avoid any modification of the data stored here:
        static void setprecomputed(std::shared_ptr<operation> op, std::vector<std::vector<densematrix>> val);
        static void setprecomputedfft(std::shared_ptr<operation> op, densematrix val);
        
        // If set to true the individual right handside contribution addresses and values are stored in 'rhsterms' when generating a formulation:
        static bool keeptrackofrhsassembly;
        // Every row in a given (int)densematrix corresponds to a shape function and every column to a given mesh element.
        // Do not forget to clear 'rhsterms' when you don't want to keep track anymore!
        static std::vector<std::pair<intdensematrix, densematrix>> rhsterms;
        
        // This stores the vec containing a solution x, its time derivative dtx and its second time derivative dtdtx
        // respectively at index 0, 1 and 2. If xdtxdtdtx[i] is an empty vector then that solution is not available.
        static std::vector<std::vector<vec>> xdtxdtdtx;
        
        
        
        // Store all !HIERARCHICAL! form function polynomials (can always be reused) and evaluated values.
        // 'formfuncpolys[i].first' gives the ith form function type name.
        // 'formfuncpolys[i].second' gives a vector detailed below.
        // 'formfuncpolys[i].second[elemtypenum][interpolorder][0]' gives the polynomials.
        static std::vector<std::pair< std::string, std::vector<std::vector< std::vector<hierarchicalformfunctioncontainer> >> >> formfuncpolys;

        // This function returns the requested form function values and reuses any already computed value if 'isreuseallowed' is true.
        // In case 'isreuseallowed' is false a pointer to the evaluated form function polynomial storage is returned for speed reasons.
        // When multiple calls follow each other and 'isreuseallowed' is false the latter storage might be modified!
        static hierarchicalformfunctioncontainer* gethff(std::string fftypename, int elementtypenumber, int interpolorder, std::vector<double> evaluationcoordinates);
        // Keep the polynomials but reset the values:
        static void resethff(void);
        
        
        // Store element split definitions. splitdefinition[elementtypenumber][splitidentifier].
        static std::vector< std::vector< std::vector<std::vector<int>> > > splitdefinition;
        // Return true if available and false otherwise.
        static bool getsplitdefinition(std::vector<std::vector<int>>& splitdef, int elementtypenumber, int splitnum, std::vector<int>& edgenumbers);
        static void setsplitdefinition(std::vector<std::vector<int>>& splitdef, int elementtypenumber, int splitnum, std::vector<int>& edgenumbers);
        
};

#endif
