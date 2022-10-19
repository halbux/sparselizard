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
#include "sl.h"
#include "universe.h"
#include "memory.h"
#include "densemat.h"
#include "indexmat.h"
#include "rawvec.h"
#include "rawmat.h"
#include "integration.h"
#include "port.h"
#include "portrelation.h"

class integration;
class contribution;
class port;
class portrelation;

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
        std::shared_ptr<dofmanager> mydofmanager = NULL;
        
        // The port relations:
        std::vector<std::shared_ptr<portrelation>> myportrelations = {};
        
        // mycontributions[m][i][j] gives the jth contribution of block number i for:
        // - the right handside if     m = 0
        // - the stiffness matrix if   m = 1
        // - the damping matrix if     m = 2
        // - the mass matrix if        m = 3
        std::vector< std::vector<std::vector<contribution>> > mycontributions = {{}, {}, {}, {}};
        
        // Always call this generate from the public generate functions:
        void generate(int m, int contributionnumber);
        
    public:
        
        formulation(void);
        
        // Add a port relation:
        formulation& operator+=(expression expr);
        
        // The following adds the contribution defined in the integration object.
        formulation& operator+=(integration integrationobject);
        formulation& operator+=(std::vector<integration> integrationobject);
        
        int countdofs(void);
        long long int allcountdofs(void);
        
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
        
        void generatein(int rhskcm, std::vector<int> contributionnumbers);
        void generate(std::vector<int> contributionnumbers);
        void generate(int contributionnumber);
        
        // Compute the no-port term value for every port relation:
        densemat getportrelationrhs(void);
        std::tuple<indexmat, indexmat, densemat> getportrelations(int KCM);
        
        std::shared_ptr<dofmanager> getdofmanager(void) { return mydofmanager; };
        
        // Get the assembled matrices or get the right hanside vector.
        // Choose to discard or not all values after getting the vector/matrix.
        
        // b() is an alias for rhs() and A() for K():
        vec b(bool keepvector = false, bool dirichletandportupdate = true);
        mat A(bool keepfragments = false);
        
        vec rhs(bool keepvector = false, bool dirichletandportupdate = true);
        mat K(bool keepfragments = false);
        mat C(bool keepfragments = false);
        mat M(bool keepfragments = false);
        // KCM set to 0 gives K, 1 gives C and 2 gives M.
        mat getmatrix(int KCM, bool keepfragments = false, std::vector<indexmat> additionalconstraints = {});
        
        
        // Generate, solve and save to fields:
        void solve(std::string soltype = "lu", bool diagscaling = false, std::vector<int> blockstoconsider = {-1});

        // DDM resolution with Dirichlet / mixed interface conditions. The initial solution is taken from the fields state. The relative residual history is returned.
        std::vector<double> allsolve(double relrestol, int maxnumit, std::string soltype = "lu", int verbosity = 1);
        std::vector<double> allsolve(std::vector<int> formulterms, std::vector<std::vector<int>> physicalterms, std::vector<std::vector<int>> artificialterms, double relrestol, int maxnumit, std::string soltype = "lu", int verbosity = 1);

};



#endif
