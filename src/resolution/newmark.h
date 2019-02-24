// sparselizard - Copyright (C) 2017- A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

// This object implements the Newmark method to solve in time the problem
//
// M*dtdtx + C*dtx + K*x = b 
//
// be it linear or nonlinear.

#ifndef NEWMARK_H
#define NEWMARK_H

#include <iostream>
#include <vector>
#include "vec.h"
#include "universe.h"
#include "mathop.h"
#include "formulation.h"

class newmark
{
    private:
        
        formulation myformulation;
        
        // The two parameters for Newmark (0.25 - 0.5 gives an unconditionally stable algorithm):
        double beta = 0.25;
        double gamma = 0.5;
        
        // The convergence tolerance for the fixed-point nonlinear iteration:
        double tol = 1e-3;
        
        // Set 'isconstant[i]' to true and the corresponding matrix/vector is 
        // supposed constant in time and will only be generated once then reused.
        //
        // - i = 0 corresponds to the rhs vector
        // - i = 1 corresponds to the K matrix
        // - i = 2 corresponds to the C matrix
        // - i = 3 corresponds to the M matrix
        //
        // Note: even if the rhs vector can be reused the Dirichlet
        // constraints will nevertheless be recomputed at each time step.
        //
        std::vector<bool> isconstant = {false, false, false, false};
        
        // The displacement u, speed v and acceleration a at the current time step:
        vec u, v, a;
        
        
        std::vector<vec> run(bool islinear, double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity);
        
    public:

        newmark(formulation formul, vec initdisplacement, vec initspeed, vec initacceleration, std::vector<bool> isrhskcmconstant = {false, false, false, false}, double b = 0.25, double g = 0.5);
        
        // Set the tolerance for the inner nonlinear fixed-point iteration:
        void settolerance(double newtol) { tol = newtol; };
        
        std::vector<vec> getcurrentsolution(void) { return {u, v, a}; };
        
        // Solve from 'starttime' to 'endtime' with constant time steps of 'timestep' 
        // seconds. One solution every 'outputeverynthtimestep' time steps is output.
        std::vector<vec> runlinear(double starttime, double timestep, double endtime, int outputeverynthtimestep = 1, int verbosity = 1);
        // Set 'maxnumnlit' to <= 0 for an unlimited number of nonlinear iterations.
        std::vector<vec> runnonlinear(double starttime, double timestep, double endtime, int maxnumnlit = -1, int outputeverynthtimestep = 1, int verbosity = 1);
        
};

#endif
