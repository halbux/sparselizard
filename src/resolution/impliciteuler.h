// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE.txt file for license information. Please report all
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
        
        // Vector x and dt(x) at the current time step:
        vec x, dtx;
        
        
        std::vector<vec> run(bool islinear, double starttime, double timestep, double endtime, int outputeverynthtimestep, int verbosity);
        
    public:

        impliciteuler(formulation formul, vec xinit, vec dtxinit, std::vector<bool> isrhskcconstant = {false, false, false});
        
        // Set the tolerance for the inner nonlinear fixed-point iteration:
        void settolerance(double newtol) { tol = newtol; };
        
        std::vector<vec> getcurrentsolution(void) { return {x, dtx}; };
        
        // Solve from 'starttime' to 'endtime' with constant time steps of 'timestep' 
        // seconds. One solution every 'outputeverynthtimestep' time steps is output.
        std::vector<vec> runlinear(double starttime, double timestep, double endtime, int outputeverynthtimestep = 1, int verbosity = 1);
        std::vector<vec> runnonlinear(double starttime, double timestep, double endtime, int outputeverynthtimestep = 1, int verbosity = 1);
        
};

#endif
