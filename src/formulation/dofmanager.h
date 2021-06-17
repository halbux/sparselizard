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
#include <unordered_map>
#include "selector.h"
#include "rawport.h"

class rawfield;
class rawport;

class dofmanager
{
    private:
        
        // The number of dofs managed:
        int numberofdofs = 0;
        
        // All the fields in the structure:
        std::vector<std::shared_ptr<rawfield>> myfields = {};
        // The currently selected field (if any) is myfields[selectedfieldnumber]:
        int selectedfieldnumber = -1;
        
        // Order of each field on each disjoint region:
        std::vector<std::vector<int>> myfieldorders = {};
        
        // Map every added port to its dof index:
        std::unordered_map<rawport*, int> myrawportmap;
        
        // 'primalondisjreg[selectedfieldnumber][disjreg]' gives
        // the primal rawport on the disjoint region (NULL if none):
        std::vector<std::vector< std::shared_ptr<rawport> >> primalondisjreg = {};
        
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
        
        bool isitmanaged = true;
        
        
        int mymeshnumber = 0;
        
        // Track the calls to 'addtostructure'.
        std::vector<std::shared_ptr<rawport>> myportstructuretracker = {};
        std::vector<std::pair<std::shared_ptr<rawfield>, int>> mystructuretracker = {};
        
        // Synchronize with the hp-adapted mesh:
        void synchronize(void);
        // To avoid infinite recursive calls:
        bool issynchronizing = false;
        
    
        // Actual function to add to the structure.
        void addtostructure(std::shared_ptr<rawfield> fieldtoadd, std::vector<int> selecteddisjointregions);
        
    public:
        
        dofmanager(void);
        // Unmanaged structure:
        dofmanager(int numdofs);
        
        bool ismanaged(void) { return isitmanaged; };
        
        void donotsynchronize(void);
        
        // Add a rawport to the structure:
        void addtostructure(std::shared_ptr<rawport> porttoadd);
        
        // 'addtostructure' defines dofs for a field on the disjoint 
        // regions. Only fields with a single component are accepted.
        void addtostructure(std::shared_ptr<rawfield> fieldtoadd, int physicalregionnumber);
        
        // Always select the field before accessing the dof structure.
        void selectfield(std::shared_ptr<rawfield> selectedfield);
        
        // Get all disjoint regions on which the selected field has dofs:
        std::vector<int> getdisjointregionsofselectedfield(void);
        
        int getrangebegin(int disjreg, int formfunc);
        int getrangeend(int disjreg, int formfunc);
        
        // Get the port dof index:
        int getaddress(rawport* prt);
        
        // Get the pointer and dof index of every port defined in this object:
        void getportsinds(std::vector<rawport*>& rps, intdensematrix& inds);
        
        // Return the primal and dual addresses for all associated ports: 
        std::pair<intdensematrix, intdensematrix> findassociatedports(void);
        
        bool isdefined(int disjreg, int formfunc);
        
        bool isported(int disjreg);
        
        // For all types of constraints:
        std::vector<bool> isconstrained(void);
        intdensematrix getconstrainedindexes(void);
        
        int countdisjregconstraineddofs(void);
        intdensematrix getdisjregconstrainedindexes(void);
        
        int countgaugeddofs(void);
        intdensematrix getgaugedindexes(void);
        
        // Get the conditionally constrained adresses as well as the constraint values:
        std::pair<intdensematrix, densematrix> getconditionalconstraintdata(void);
        
        std::shared_ptr<rawfield> getselectedfield(void);
        std::vector<std::shared_ptr<rawfield>> getfields(void);
        // The replacing field must have an identical type and order vector as the replaced one:
        void replaceselectedfield(std::shared_ptr<rawfield> rf);
        
        std::vector<int> getselectedfieldorders(void);
        
        // Count the total number of ports (primal, dual and not-associated):
        int countports(void);
        int countassociatedprimalports(void);
        
        int countdofs(void);
        long long int allcountdofs(void);
        int countformfunctions(int disjointregion);
        
        // Return {sendnewconstrainedinds, recvnewconstrainedinds, sendunconstrainedinds, recvunconstrainedinds} where
        //
        // - sendnewconstrainedinds[n] are the indexes of all interface dofs constrained on this rank but not constrained on the neighbour
        // - recvnewconstrainedinds[n] are the indexes of all interface dofs not constrained on this rank but constrained on the neighbour
        // - sendunconstrainedinds[n] are all unconstrained indexes in senddofinds (same ordering)
        // - recvunconstrainedinds[n] are all unconstrained indexes in recvdofinds (same ordering)
        //
        std::vector<std::vector<intdensematrix>> discovernewconstraints(std::vector<int> neighbours, std::vector<intdensematrix> senddofinds, std::vector<intdensematrix> recvdofinds);
        
        void print(void);
        
        // 'getaddresses' is required in the matrix generation step.
        // It returns an intdensematrix representing a numberofformfunctions
        // by elementlist.size() matrix (row-major). The matrix gives 
        // the addresses in the formulation matrix at which the dofs of field 
        // 'inputfield' defined on the elements in elementlist can be found. 
        // Address -1 is used for field dofs not in 'fieldphysreg'.
        //
        intdensematrix getaddresses(std::shared_ptr<rawfield> inputfield, int fieldinterpolationorder, int elementtypenumber, std::vector<int> &elementlist, int fieldphysreg);
                                                        
};

#endif
