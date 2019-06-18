// sparselizard - Copyright (C) 2017- A. Halbach
//
// Contributions kindly provided by R. Haouari.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

// This object implements the generalized alpha method to solve in time the problem
//
// M*dtdtx + C*dtx + K*x = b 
//
// be it linear or nonlinear.

#ifndef GENALPHA_H
#define GENALPHA_H

#include <iostream>
#include <vector>
#include "vec.h"
#include "universe.h"
#include "mathop.h"
#include "formulation.h"

class genalpha
{
    private:
        
        formulation myformulation;
        
        // The four parameters for generalized alpha (set to unconditionally stable Newmark by default):
        double beta = 0.25, gamma = 0.5, alphaf = 0.0, alpham = 0.0;
        
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
        
        // Formulations to solve before/after 'myformulation' is solved:
        std::vector<formulation> tosolvebefore = {};
        std::vector<formulation> tosolveafter = {};
        
        // The displacement u, speed v and acceleration a at the current time step:
        vec u, v, a;
        
        
        std::vector<std::vector<vec>> run(bool islinear, double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity);
        
    public:
    
        genalpha(formulation formul, vec initdisplacement, vec initspeed, vec initacceleration, std::vector<bool> isrhskcmconstant = {false, false, false, false});
        
        // Manually specify all four parameters:
        void setparameter(double b, double g, double af, double am) { beta = b; gamma = g; alphaf = af; alpham = am; };
        // Specify a high-frequency dissipation and let the four parameters be optimally deduced:
        void setparameter(double rinf);
        
        // Set the tolerance for the inner nonlinear fixed-point iteration:
        void settolerance(double newtol) { tol = newtol; };
        
        std::vector<vec> getcurrentsolution(void) { return {u, v, a}; };
        
        // Define a list of formulations to solve at the beginning/end of the nonlinear loop:
        void presolve(std::vector<formulation> formuls);
        void postsolve(std::vector<formulation> formuls);
        
        // Solve from 'starttime' to 'endtime' with constant time steps of 'timestep' 
        // seconds. output[0] gives the u time-solutions while output[1] gives v and output[2] gives a. 
        // One solution every 'outputeverynthtimestep' time steps is output.
        std::vector<std::vector<vec>> runlinear(double starttime, double timestep, double endtime, int outputeverynthtimestep = 1, int verbosity = 1);
        // Set 'maxnumnlit' to <= 0 for an unlimited number of nonlinear iterations.
        std::vector<std::vector<vec>> runnonlinear(double starttime, double timestep, double endtime, int maxnumnlit = -1, int outputeverynthtimestep = 1, int verbosity = 1);
        
};

#endif
