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

void impliciteuler::presolve(std::vector<formulation> formuls) { tosolvebefore = formuls; }
void impliciteuler::postsolve(std::vector<formulation> formuls) { tosolveafter = formuls; }

std::vector<std::vector<vec>> impliciteuler::runlinear(double starttime, double timestep, double endtime, int outputeverynthtimestep, int verbosity)
{
    return run(true, starttime, timestep, endtime, -1, outputeverynthtimestep, verbosity);
}

std::vector<std::vector<vec>> impliciteuler::runnonlinear(double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity)
{
    return run(false, starttime, timestep, endtime, maxnumnlit, outputeverynthtimestep, verbosity);
}

std::vector<std::vector<vec>> impliciteuler::run(bool islinear, double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity)
{
    // Solve end time rounding issues:
    endtime += endtime*1e-12;
    
    if (starttime > endtime)
        return {};
    
    if (outputeverynthtimestep <= 0)
        outputeverynthtimestep = 1;
    
    // Make all time derivatives available in the universe:
    universe::xdtxdtdtx = {{x},{dtx},{}};
        
    // Get all fields in the formulation:
    std::vector<std::shared_ptr<rawfield>> allfields = myformulation.getdofmanager()->getfields();
    // Set all fields in the formulation to the initial solution:
    for (int i = 0; i < allfields.size(); i++)
        field(allfields[i]).setdata(-1, x);
    
    vec rhs; mat K, C, leftmat;
    
    // Count the number of time steps to step through and the number of vectors to output:
    int numtimesteps = 0; int outputsize = 0;
    for (double t = starttime; t <= endtime; t = t + timestep)
    {
        if (numtimesteps%outputeverynthtimestep == 0)
            outputsize++;
        numtimesteps++;
    }
    
    
    // Start the implicit Euler iteration:
    std::cout << "Implicit Euler for " << numtimesteps << " timesteps in range " << starttime << " to " << endtime << " sec:" << std::endl;
    std::vector<std::vector<vec>> output(2, std::vector<vec>(outputsize));
    output[0][0] = x; output[1][0] = dtx;
    
    // We already have everything for time step 0 so we start at 1:
    int timestepindex = 1;
    for (double t = starttime + timestep; t <= endtime; t = t + timestep)
    {        
        std::cout << timestepindex << "@" << t << "sec" << std::flush;

        mathop::settime(t);
        
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
            if (isconstant[1] == false || timestepindex == 1)
            {
                myformulation.generatestiffnessmatrix();
                K = myformulation.K(false, false);
            }
            if (isconstant[2] == false || timestepindex == 1)
            {
                myformulation.generatedampingmatrix();
                C = myformulation.C(false, true);
            }
            if (isconstant[0] == false || timestepindex == 1)
            {
                myformulation.generaterhs();
                rhs = myformulation.rhs();
            }
            else
                rhs.updateconstraints();
            
            // Reuse matrices when possible (including the LU decomposition):
            if (isconstant[1] == false || isconstant[2] == false || timestepindex == 1)
            {
                leftmat = C + timestep*K;
                leftmat.reuselu();
            }
            
            // Update the solution xnext.
            xnext = relaxationfactor * mathop::solve(leftmat, C*x+timestep*rhs) + (1.0-relaxationfactor)*xnext;
            
            dtxnext = 1.0/timestep*(xnext-x);
            
            // Update all fields in the formulation:
            for (int i = 0; i < allfields.size(); i++)
                field(allfields[i]).setdata(-1, xnext);
            
            relchange = (xnext-xtolcalc).norm()/xnext.norm();
            
            if (islinear == false && verbosity > 0)
                std::cout << " " << relchange << std::flush;

            nlit++; 
            
            
            // Solve all formulations that must be solved at the end of the nonlinear loop:
            mathop::solve(tosolveafter);
            
            
            if (islinear)
                break;
        }
        
        x = xnext;
        dtx = dtxnext;
        
        if (islinear == false)
            std::cout << " (" << nlit << "NL it)" << std::flush;
        if (timestepindex < numtimesteps-1)
        std::cout << " -> " << std::flush;
        
        // Only one every 'outputeverynthtimestep' solutions is output:
        if (timestepindex%outputeverynthtimestep == 0)
        {
            output[0][timestepindex/outputeverynthtimestep] = x;
            output[1][timestepindex/outputeverynthtimestep] = dtx;
        }
        timestepindex++;
    }
    std::cout << std::endl;
    
    // Remove all time derivatives from the universe:
    universe::xdtxdtdtx = {{},{},{}};
    
    return output;
}


