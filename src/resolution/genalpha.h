// sparselizard - Copyright (C) see copyright file.
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
        
        int myverbosity = 1;
        
        formulation myformulation;
        
        // The four parameters for generalized alpha (set to unconditionally stable Newmark by default):
        double beta = 0.25, gamma = 0.5, alphaf = 0.0, alpham = 0.0;
        
        // The convergence tolerance for the fixed-point nonlinear iteration:
        double nltol = 1e-3;
        
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
        
        // Current timestep:
        double dt = -1;
        // All time values stepped-through:
        std::vector<double> mytimes = {};
        
        // Time-adaptivity settings:
        double mindt = -1, maxdt = -1, tatol = -1, rfact = -1, cfact = -1, cthres = -1;
        
        // The speed v and acceleration a at the current time step:
        vec v, a;
        
        // Objects required at every timestep (possibly reused):
        vec rhs; mat K, C, M, leftmat, matu, matv, mata;
        // Parameters for which these objects are defined:
        double defbeta = -1, defgamma = -1, defalphaf = -1, defalpham = -1, defdt = -1;
        
        int run(bool islinear, double timestep, int maxnumnlit);
        
    public:
    
        genalpha(formulation formul, vec initspeed, vec initacceleration, int verbosity = 3, std::vector<bool> isrhskcmconstant = {false, false, false, false});
    
        void setverbosity(int verbosity) { myverbosity = verbosity; };
        
        // Manually specify all four parameters:
        void setparameter(double b, double g, double af, double am) { beta = b; gamma = g; alphaf = af; alpham = am; };
        // Specify a high-frequency dissipation and let the four parameters be optimally deduced:
        void setparameter(double rinf);
        
        // Set the tolerance for the inner nonlinear fixed-point iteration:
        void settolerance(double tol) { nltol = tol; };
        
        std::vector<vec> gettimederivative(void) { return {v, a}; };
        void settimederivative(std::vector<vec> sol);
        
        void settimestep(double timestep) { dt = timestep; };
        double gettimestep(void) { return dt; };
        
        // Count the number of timesteps computed:
        int count(void) { return mytimes.size(); };
        std::vector<double> gettimes(void) { return mytimes; };
        
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
