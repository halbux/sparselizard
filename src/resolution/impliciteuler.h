// sparselizard - Copyright (C) 2017- A. Halbach
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
        double tol = 1e-3;
        
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
        
        // Vector x and dt(x) at the current time step:
        vec x, dtx;
        
        
        std::vector<std::vector<vec>> run(bool islinear, double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity);
        
    public:

        impliciteuler(formulation formul, vec xinit, vec dtxinit, std::vector<bool> isrhskcconstant = {false, false, false});
        
        // Set the tolerance for the inner nonlinear fixed-point iteration:
        void settolerance(double newtol) { tol = newtol; };
        
        // Set the relaxation factor for the inner nonlinear fixed-point iteration:
        void setrelaxationfactor(double relaxfact) { relaxationfactor = relaxfact; };
        
        std::vector<vec> getcurrentsolution(void) { return {x, dtx}; };
        
        // Define a list of formulations to solve at the beginning/end of the nonlinear loop:
        void presolve(std::vector<formulation> formuls);
        void postsolve(std::vector<formulation> formuls);
        
        // Solve from 'starttime' to 'endtime' with constant time steps of 'timestep' 
        // seconds. output[0] gives the x time-solutions while output[1] gives dt(x). 
        // One solution every 'outputeverynthtimestep' time steps is output.
        std::vector<std::vector<vec>> runlinear(double starttime, double timestep, double endtime, int outputeverynthtimestep = 1, int verbosity = 1);
        // Set 'maxnumnlit' to <= 0 for an unlimited number of nonlinear iterations.
        std::vector<std::vector<vec>> runnonlinear(double starttime, double timestep, double endtime, int maxnumnlit = -1, int outputeverynthtimestep = 1, int verbosity = 1);
        
};

#endif
