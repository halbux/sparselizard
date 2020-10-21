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
        
        // Vector dt(x) at the current time step:
        vec dtx;
        
        // Objects required at every timestep (possibly reused):
        vec rhs; mat K, C, leftmat;
        // Parameters for which these objects are defined:
        double defdt = -1;
        
        int run(bool islinear, double timestep, int maxnumnlit, int verbosity, bool autoadvancetime);
        
    public:

        impliciteuler(formulation formul, vec dtxinit, std::vector<bool> isrhskcconstant = {false, false, false});
        
        // Set the tolerance for the inner nonlinear fixed-point iteration:
        void settolerance(double newtol) { nltol = newtol; };
        
        // Set the relaxation factor for the inner nonlinear fixed-point iteration:
        void setrelaxationfactor(double relaxfact) { relaxationfactor = relaxfact; };
        
        vec gettimederivative(void) { return dtx; };
        void settimederivative(vec sol);
        
        // Get the timestep:
        double gettimestep(void) { return dt; };
        
        // Define a list of formulations to solve at the beginning/end of the nonlinear loop:
        void presolve(std::vector<formulation> formuls);
        void postsolve(std::vector<formulation> formuls);
        
        // Advance the solution by the provided timestep.
        void runlinear(double timestep, int verbosity = 1, bool autoadvancetime = true);
        // Set 'maxnumnlit' to <= 0 for an unlimited number of nonlinear iterations.
        int runnonlinear(double timestep, int maxnumnlit = -1, int verbosity = 2, bool autoadvancetime = true);
        
};

#endif
