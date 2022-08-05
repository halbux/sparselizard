// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

// This code calls the PETSc library. See https://www.mcs.anl.gov/petsc/ for more information.

#ifndef RAWVEC_H
#define RAWVEC_H

#include <iostream>
#include <string>
#include "dofmanager.h"
#include "indexmat.h"
#include "densemat.h"
#include "vectorfieldselect.h"
#include "rawfield.h"
#include "memory.h"
#include "petsc.h"
#include "petscvec.h"
#include "sl.h"
#include "ptracker.h"
#include "rawmesh.h"

class vectorfieldselect;
class dofmanager;
class rawfield;
class rawmesh;

class rawvec : public std::enable_shared_from_this<rawvec>
{
    private:

        Vec myvec = PETSC_NULL;
        std::shared_ptr<dofmanager> mydofmanager = NULL;
        

        std::shared_ptr<ptracker> myptracker = NULL;
        std::vector<dofmanager> mycurrentstructure = {};
        
        // Synchronize with the hp-adapted mesh:
        void synchronize(void);
        // To avoid infinite recursive calls:
        bool issynchronizing = false;
        // Allow/forbid value syncing:
        bool isvaluesynchronizingallowed = true;
        
        
        // Mesh on which this object is based:
        std::shared_ptr<rawmesh> myrawmesh = NULL;
    
    public:
            
        rawvec(std::shared_ptr<dofmanager> dofmngr);
        rawvec(std::shared_ptr<dofmanager> dofmngr, Vec input);
        
        ~rawvec(void);
        
        void allowvaluesynchronizing(bool allowit);
        
        int size(void);
        
        // Update the indexes that correspond to constrained values of a rawfield:
        void updatedisjregconstraints(std::shared_ptr<rawfield> constrainedfield);
        
        // Negative addresses are ignored. 'op' can be 'add' or 'set'. 
        void setvalues(indexmat addresses, densemat valsmat, std::string op = "set");
        // Set value at a single index:
        void setvalue(int address, double value, std::string op = "set");
        
        densemat getvalues(indexmat addresses);
        // Get value at a single index:
        double getvalue(int address);
        
        void setvalues(std::shared_ptr<rawfield> selectedfield, int disjointregionnumber, int formfunctionindex, densemat vals, std::string op);
        densemat getvalues(std::shared_ptr<rawfield> selectedfield, int disjointregionnumber, int formfunctionindex);
        
        void setvaluestoports(void);
        void setvaluesfromports(void);
        
        // Write and load raw vec data:
        void write(std::string filename);
        void load(std::string filename);
        
        void print(void);
        
        std::shared_ptr<dofmanager> getdofmanager(void);
        Vec getpetsc(void);
        
        std::shared_ptr<rawvec> getpointer(void);
        std::shared_ptr<rawmesh> getrawmesh(void);
        std::shared_ptr<ptracker> getptracker(void);
        
        // Transfer data between the 'inputvec' vector and this vector:
        void setdata(std::shared_ptr<rawvec> inputvec, int disjreg, std::shared_ptr<rawfield> inputfield);
        
};

#endif
