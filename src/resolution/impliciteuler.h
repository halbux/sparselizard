// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

// This object implements the implicit Euler (also called backward Euler) method to solve in time the problem
//
// C*dtx + K*x = b 
//
// be it linear or nonlinear.

#ifndef IMPLICITEULER_H
#define IMPLICITEULER_H

#include <iostream>
#include <vector>
#include "vec.h"
#include "universe.h"
#include "mathop.h"
#include "formulation.h"

class impliciteuler
{
    private:
        
        int myverbosity = 1;
        
        formulation myformulation;
        
        // The convergence tolerance for the fixed-point nonlinear iteration:
        double nltol = 1e-3;
        
        // The relaxation factor for the nonlinear iteration:
        double relaxationfactor = 1.0;
        
        // Set 'isconstant[i]' to true and the corresponding matrix/vector is 
        // supposed constant in time and will only be generated once then reused.
        //
        // - i = 0 corresponds to the rhs vector
        // - i = 1 corresponds to the K matrix
        // - i = 2 corresponds to the C matrix
        //
        // Note: even if the rhs vector can be reused the Dirichlet
        // constraints will nevertheless be recomputed at each time step.
        //
        std::vector<bool> isconstant = {false, false, false};
        
        // Formulations to solve before/after 'myformulation' is solved:
        std::vector<formulation> tosolvebefore = {};
        std::vector<formulation> tosolveafter = {};
        
        // Current timestep:
        double dt = -1;
        // Time-adaptivity settings:
        double mindt = -1, maxdt = -1, tatol = -1, rfact = -1, cfact = -1, cthres = -1;
        
        // Vector dt(x) at the current time step:
        vec dtx;
        
        // Objects required at every timestep (possibly reused):
        vec rhs; mat K, C, leftmat;
        // Parameters for which these objects are defined:
        double defdt = -1;
        
        int run(bool islinear, double timestep, int maxnumnlit);
        
    public:

        impliciteuler(formulation formul, vec dtxinit, int verbosity = 3, std::vector<bool> isrhskcconstant = {false, false, false});
        
        void setverbosity(int verbosity) { myverbosity = verbosity; };
        
        // Set the tolerance for the inner nonlinear fixed-point iteration:
        void settolerance(double tol) { nltol = tol; };
        
        // Set the relaxation factor for the inner nonlinear fixed-point iteration:
        void setrelaxationfactor(double relaxfact) { relaxationfactor = relaxfact; };
        
        vec gettimederivative(void) { return dtx; };
        void settimederivative(vec sol);
        
        void settimestep(double timestep) { dt = timestep; };
        double gettimestep(void) { return dt; };

        // Set the time-adaptivity settings:
        void setadaptivity(double tol, double mints, double maxts, double reffact = 0.5, double coarfact = 2.0, double coarthres = 0.5);
        
        // Define a list of formulations to solve at the beginning/end of the nonlinear loop:
        void presolve(std::vector<formulation> formuls);
        void postsolve(std::vector<formulation> formuls);
        
        // Advance the solution by the provided timestep for a linear/nonlinear problem.
        void next(double timestep);
        int next(double timestep, int maxnumnlit);
        
};

#endif
