#include "impliciteuler.h"

impliciteuler::impliciteuler(formulation formul, vec xinit, vec dtxinit, std::vector<bool> isrhskcconstant)
{
    myformulation = formul;
    
    x = xinit;
    dtx = dtxinit;
    
    if (isrhskcconstant.size() == 0)
        isconstant = {false,false,false};
    else
        isconstant = isrhskcconstant;
    if (isconstant.size() != 3)
    {
        std::cout << "Error in 'impliciteuler' object: expected a length 3 or empty vector as fourth argument" << std::endl;
        abort();  
    }
}

void impliciteuler::setsolution(std::vector<vec> sol)
{
    if (sol.size() != 2)
    {
        std::cout << "Error in 'impliciteuler' object: expected a vector of length two to set the solution" << std::endl;
        abort();  
    }
    x = sol[0]; dtx = sol[1];
}

void impliciteuler::presolve(std::vector<formulation> formuls) { tosolvebefore = formuls; }
void impliciteuler::postsolve(std::vector<formulation> formuls) { tosolveafter = formuls; }

void impliciteuler::runlinear(double timestep, int verbosity, bool autoadvancetime)
{
    run(true, timestep, -1, verbosity, autoadvancetime);
}

int impliciteuler::runnonlinear(double timestep, int maxnumnlit, int verbosity, bool autoadvancetime)
{
    return run(false, timestep, maxnumnlit, verbosity, autoadvancetime);
}

int impliciteuler::run(bool islinear, double timestep, int maxnumnlit, int verbosity, bool autoadvancetime)
{
    double inittime = universe::currenttimestep;

    dt = timestep;
    // Update and print the time:
    universe::currenttimestep += dt;
    char spacer = ':';
    if (islinear || verbosity < 2)
        spacer = ' ';
    if (verbosity > 0)
        std::cout << "@" << universe::currenttimestep << "s" << spacer << std::flush;
    
    // Make all time derivatives available in the universe:
    universe::xdtxdtdtx = {{x},{dtx},{}};
        
    // Set all fields in the formulation to the initial solution:
    mathop::setdata(x);
    
    // Nonlinear loop:
    double relchange = 1; int nlit = 0;
    vec xnext = x, dtxnext = dtx;
    while (relchange > tol && (maxnumnlit <= 0 || nlit < maxnumnlit))
    {
        // Solve all formulations that must be solved at the beginning of the nonlinear loop:
        mathop::solve(tosolvebefore);

        // Make all time derivatives available in the universe:
        universe::xdtxdtdtx = {{xnext},{dtxnext},{}};
        
        vec xtolcalc = xnext;
        
        // Reassemble only the non-constant matrices:
        bool isfirstcall = not(K.isdefined());
        if (isconstant[1] == false || isfirstcall)
        {
            myformulation.generatestiffnessmatrix();
            K = myformulation.K(false, false);
        }
        if (isconstant[2] == false || isfirstcall)
        {
            myformulation.generatedampingmatrix();
            C = myformulation.C(false, true);
        }
        if (isconstant[0] == false || isfirstcall)
        {
            myformulation.generaterhs();
            rhs = myformulation.rhs();
        }
        else
            rhs.updateconstraints();
        
        // Reuse matrices when possible (including the factorization):
        if (isconstant[1] == false || isconstant[2] == false || isfirstcall || defdt != dt)
        {
            leftmat = C + dt*K;
            leftmat.reusefactorization();
            
            defdt = dt;
        }
        
        // Update the solution xnext.
        xnext = relaxationfactor * mathop::solve(leftmat, C*x+dt*rhs) + (1.0-relaxationfactor)*xnext;
        
        dtxnext = 1.0/dt*(xnext-x);
        
        // Update all fields in the formulation:
        mathop::setdata(xnext);
        
        relchange = (xnext-xtolcalc).norm()/xnext.norm();
        
        if (islinear == false && verbosity > 1)
            std::cout << " " << relchange << std::flush;

        nlit++; 
        
        // Solve all formulations that must be solved at the end of the nonlinear loop:
        mathop::solve(tosolveafter);
        
        if (islinear)
            break;
    }
    
    x = xnext; dtx = dtxnext;
    
    if (verbosity > 1 && islinear == false)
        std::cout << " (" << nlit << "NL it) " << std::flush;
    
   if (autoadvancetime == false)
        universe::currenttimestep = inittime;
    
    return nlit;
}


