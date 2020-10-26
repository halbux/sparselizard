// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef FORMULATION_H
#define FORMULATION_H

#include <iostream>
#include <vector>
#include "integration.h"
#include "contribution.h"
#include "dofmanager.h"
#include "vec.h"
#include "mat.h"
#include "expression.h"
#include "mathop.h"
#include "universe.h"
#include "memory.h"
#include "densematrix.h"
#include "intdensematrix.h"
#include "rawvec.h"
#include "rawmat.h"
#include "integration.h"

class integration;
class contribution;

class formulation
{
    private:
        
        bool isstructurelocked = false;
        
        // myvec is the right handside vector rhs.
        std::shared_ptr<rawvec> myvec = NULL;
        // - mymat[0] is the stiffness matrix K
        // - mymat[1] is the damping matrix C
        // - mymat[2] is the mass matrix M
        std::vector<std::shared_ptr<rawmat>> mymat = {NULL, NULL, NULL};
        
        // The link between the dof number and its row and column in the matrix:
        std::shared_ptr<dofmanager> mydofmanager;
        
        // mycontributions[m][i][j] gives the jth contribution of block number i for:
        // - the right handside if     m = 0
        // - the stiffness matrix if   m = 1
        // - the damping matrix if     m = 2
        // - the mass matrix if        m = 3
        std::vector< std::vector<std::vector<contribution>> > mycontributions = {{}, {}, {}, {}};
        
        // Always call this generate from the public generate functions:
        void generate(int m, int contributionnumber);
        
    public:
        
        // Has this formulation been called to compute a constraint?
        bool isconstraintcomputation = false;
        
        
        formulation(void);
        
        // The following adds the contribution defined in the integration object.
        void operator+=(integration integrationobject);
        void operator+=(std::vector<integration> integrationobject);
        
        int countdofs(void);
        
        bool isstiffnessmatrixdefined(void);
        bool isdampingmatrixdefined(void);
        bool ismassmatrixdefined(void);
        
        // Generate all blocks:
        void generate(void);
        // Generate all contributions of K, C, M or rhs:
        void generatestiffnessmatrix(void);
        void generatedampingmatrix(void);
        void generatemassmatrix(void);
        void generaterhs(void);
        
        void generate(std::vector<int> contributionnumbers);
        void generate(int contributionnumber);
        
        std::shared_ptr<dofmanager> getdofmanager(void) { return mydofmanager; };
        
        // Get the assembled matrices or get the right hanside vector.
        // Choose to discard or not all values after getting the vector/matrix.
        // Choose to add or not the diagonal ones for the Dirichlet constraints.
        
        // b() is an alias for rhs() and A() for K():
        vec b(bool keepvector = false);
        mat A(bool keepfragments = false, bool skipdiagonalones = false);
        
        vec rhs(bool keepvector = false);
        mat K(bool keepfragments = false, bool skipdiagonalones = false);
        mat C(bool keepfragments = false, bool skipdiagonalones = false);
        mat M(bool keepfragments = false, bool skipdiagonalones = false);
        // KCM set to 0 gives K, 1 gives C and 2 gives M.
        mat getmatrix(int KCM, bool keepfragments = false, bool skipdiagonalones = false);

};



#endif
