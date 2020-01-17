// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// The 'dofmanager' object manages the structure which defines at which 
// indexes of an assembled matrix the degrees of freedom can be found.

#ifndef DOFMANAGER_H
#define DOFMANAGER_H

#include "rawfield.h"
#include "element.h"
#include "universe.h"
#include <vector>
#include "disjointregions.h"
#include "intdensematrix.h"
#include <memory>
#include "selector.h"

class rawfield;

class dofmanager
{
    private:
        
        // The number of dofs managed:
        int numberofdofs = 0;
        
        // All the fields in the structure:
        std::vector<std::shared_ptr<rawfield>> myfields = {};
        // The currently selected field (if any) is myfields[selectedfieldnumber]:
        int selectedfieldnumber = -1;
        
        // 'rangebegin[selectedfieldnumber][12][2]' gives the index of 
        // the first row/column in the matrix at which the data for 
        //
        // - field pointed by myfields[selectedfieldnumber]
        // - vertex, edge, face or volume type-disjoint region 12
        // - the third vertex, edge, face or volume form function
        //
        // can be found. 
        // 
        // 'rangeend[selectedfieldnumber][12][2]' gives the last row.
        //
        std::vector<std::vector<std::vector< int >>> rangebegin = {};
        std::vector<std::vector<std::vector< int >>> rangeend = {};
        
        // Actual function to add to the structure. Must be private.
        void addtostructure(std::shared_ptr<rawfield> fieldtoadd, std::vector<int> selecteddisjointregions);
        
    public:
        
        dofmanager(void) {};
        dofmanager(int numdofs) { numberofdofs = numdofs; };
        
        // 'addtostructure' defines dofs for a field on the disjoint 
        // regions. Only fields with a single component are accepted.
        void addtostructure(std::shared_ptr<rawfield> fieldtoadd, int physicalregionnumber);
        // Always select the field before accessing the dof structure.
        void selectfield(std::shared_ptr<rawfield> selectedfield);
        
        // Get all disjoint regions on which the selected field has dofs:
        std::vector<int> getdisjointregionsofselectedfield(void);
        
        int getrangebegin(int disjreg, int formfunc);
        int getrangeend(int disjreg, int formfunc);
        
        bool isdefined(int disjreg, int formfunc) { return (formfunc < rangebegin[selectedfieldnumber][disjreg].size()); };
        
        int countconstraineddofs(void);
        intdensematrix getconstrainedindexes(void);
        
        int countgaugeddofs(void);
        intdensematrix getgaugedindexes(void);
        
        // Get the conditionally constrained adresses as well as the constraint values:
        std::pair<intdensematrix, densematrix> getconditionalconstraintdata(void);
        
        // Return a new dofmanager object that does not include the Dirichlet constraints.
        // 'dofrenumbering' must have a size equal to the number of dofs before the call.
        // Removed dofs are renumbered as -1.
        std::shared_ptr<dofmanager> removeconstraints(int* dofrenumbering);
        
        std::vector<std::shared_ptr<rawfield>> getfields(void) { return myfields; };
        
        int countdofs(void) { return numberofdofs; };
        
        void print(void);
        
        // 'getadresses' is required in the matrix generation step.
        // It returns an intdensematrix representing an numberofformfunctions
        // by elementlist.size() matrix (row-major). The matrix gives 
        // the adresses in the formulation matrix at which the dofs of field 
        // 'inputfield' defined on the elements in elementlist can be found. 
        //
        // The following adress tags have a special meaning:
        //
        // - adress -1 is used for constrained fields (requires 'useminusonetag' true)
        // - adress -2 is used for field dofs not in 'fieldphysreg' 
        //
        intdensematrix getadresses(std::shared_ptr<rawfield> inputfield, int fieldinterpolationorder, int elementtypenumber, std::vector<int> &elementlist, int fieldphysreg, bool useminusonetag);
                                                        
};

#endif
