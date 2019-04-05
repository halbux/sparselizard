#include "newmark.h"

newmark::newmark(formulation formul, vec initdisplacement, vec initspeed, vec initacceleration, std::vector<bool> isrhskcmconstant, double b, double g)
{
    myformulation = formul;
    
    beta = b;
    gamma = g;
    
    u = initdisplacement;
    v = initspeed;
    a = initacceleration;
    
    if (isrhskcmconstant.size() == 0)
        isconstant = {false,false,false,false};
    else
        isconstant = isrhskcmconstant;
    if (isconstant.size() != 4)
    {
        std::cout << "Error in 'newmark' object: expected a length 4 or empty vector as fourth argument" << std::endl;
        abort();  
    }
    
    // Due to a computation below beta cannot be zero:
    if (std::abs(beta) < 0.01)
    {
        std::cout << "Error in 'newmark' object: parameter beta = 0 is not supported (too close to zero should be avoided too)" << std::endl;
        abort();  
    }
}

std::vector<std::vector<vec>> newmark::runlinear(double starttime, double timestep, double endtime, int outputeverynthtimestep, int verbosity)
{
    return run(true, starttime, timestep, endtime, -1, outputeverynthtimestep, verbosity);
}

std::vector<std::vector<vec>> newmark::runnonlinear(double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity)
{
    return run(false, starttime, timestep, endtime, maxnumnlit, outputeverynthtimestep, verbosity);
}

std::vector<std::vector<vec>> newmark::run(bool islinear, double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity)
{
    // Solve end time rounding issues:
    endtime += endtime*1e-12;
    
    if (starttime > endtime)
        return {};
    
    if (outputeverynthtimestep <= 0)
        outputeverynthtimestep = 1;
        
    // Get all fields in the formulation:
    shared_ptr<dofmanager> dofmngr = myformulation.getdofmanager();
    std::vector<shared_ptr<rawfield>> allfields = dofmngr->getfields();
    // Set all fields in the formulation to the initial displacement:
    for (int i = 0; i < allfields.size(); i++)
        allfields[i]->setdata(-1, u|field(allfields[i]));
    
    // Get all indexes at which the fields are constrained:
    intdensematrix constraintindexes = myformulation.getdofmanager()->getconstrainedindexes();
    
    // Remove leftovers (if any):
    myformulation.rhs(); myformulation.K(); myformulation.C(); myformulation.M();

    vec rhs; mat K, C, M, leftmat, matu, matv, mata;
    
    // Count the number of time steps to step through and the number of vectors to output:
    int numtimesteps = 0; int outputsize = 0;
    for (double t = starttime; t <= endtime; t = t + timestep)
    {
        if (numtimesteps%outputeverynthtimestep == 0)
            outputsize++;
        numtimesteps++;
    }
    
    
    // Start the Newmark iteration:
    std::cout << "Newmark (beta " << beta << ", gamma " << gamma << ") for " << numtimesteps << " timesteps in range " << starttime << " to " << endtime << " sec:" << std::endl;
    std::vector<std::vector<vec>> output(3, std::vector<vec>(outputsize));
    output[0][0] = u; output[1][0] = v; output[2][0] = a;
    
    // We already have everything for time step 0 so we start at 1:
    int timestepindex = 1;
    for (double t = starttime + timestep; t <= endtime; t = t + timestep)
    {        
        std::cout << timestepindex << "@" << t << "sec" << std::flush;

        mathop::settime(t);
        
        // Nonlinear loop:
        double relchange = 1; int nlit = 0;
        vec unext = u, vnext = v, anext = a;
        while (relchange > tol && (maxnumnlit <= 0 || nlit < maxnumnlit))
        {
        	// Make all time derivatives available in the universe:
        	universe::xdtxdtdtx = {{unext},{vnext},{anext}};
        
            vec utolcalc = unext;
            
            // Reassemble only the non-constant matrices:
            if (isconstant[1] == false || timestepindex == 1)
            {
                myformulation.generatestiffnessmatrix();
                K = myformulation.K(false, true);
            }
            if (isconstant[2] == false || timestepindex == 1)
            {
                myformulation.generatedampingmatrix();
                C = myformulation.C(false, true);
            }
            if (isconstant[3] == false || timestepindex == 1)
            {
                myformulation.generatemassmatrix();
                M = myformulation.M(false, false);
            }
            if (isconstant[0] == false || timestepindex == 1)
            {
                myformulation.generaterhs();
                rhs = myformulation.rhs();
            }
            else
                rhs.updateconstraints();
            
            // Reuse matrices when possible (including the LU decomposition):
            if (isconstant[1] == false || isconstant[2] == false || isconstant[3] == false || timestepindex == 1)
            {
                leftmat = M + (gamma*timestep)*C + (beta*timestep*timestep)*K;
                leftmat.reuselu();
                
                matu = -K;
                matv = -(C+timestep*K);
                mata = ((gamma-1)*timestep)*C+((beta-0.5)*timestep*timestep)*K;
            }
            
            // Update the acceleration. 
            // The acceleration is imposed on the Dirichlet-constrained dofs.
            // The displacement update relation is used to make sure the 
            // acceleration constraint leads to the exact constrained 
            // displacement at the next time step:
            vec unextdirichlet(myformulation); 
            unextdirichlet.updateconstraints();
            vec anextdirichlet = 1/(beta*timestep*timestep)*( unextdirichlet-u - timestep*v - timestep*timestep*(0.5-beta)*a ); // beta cannot be zero!
            // Here are the constrained values of the next acceleration:
            densematrix anextdirichletval = anextdirichlet.getpointer()->getvalues(constraintindexes);
            // Here for the conditional constraints:
            intdensematrix condconstrainedindexes = (myformulation.getdofmanager()->getconditionalconstraintdata()).first;
            densematrix condconstranextdirichletval = anextdirichlet.getpointer()->getvalues(condconstrainedindexes);
            
            vec rightvec = matu*u + matv*v + mata*a + rhs;
            // Force the acceleration on the constrained dofs:
            rightvec.getpointer()->setvalues(constraintindexes, anextdirichletval);
            rightvec.getpointer()->setvalues(condconstrainedindexes, condconstranextdirichletval);
            
            anext = mathop::solve(leftmat, rightvec);

            // Update unext and vnext:
            unext = u + timestep*v + ((0.5-beta)*timestep*timestep)*a + (beta*timestep*timestep)*anext;
            vnext = v + (timestep*(1-gamma))*a + (gamma*timestep)*anext;
            
            // Update all fields in the formulation:
            for (int i = 0; i < allfields.size(); i++)
                allfields[i]->setdata(-1, unext|field(allfields[i]));
            
            relchange = (unext-utolcalc).norm()/unext.norm();
            
            if (islinear == false && verbosity > 0)
                std::cout << " " << relchange << std::flush;

            nlit++; 
            
            if (islinear)
                break;
        }

        u = unext; v = vnext; a = anext;
        
        if (islinear == false)
            std::cout << " (" << nlit << "NL it)" << std::flush;
        if (timestepindex < numtimesteps-1)
        std::cout << " -> " << std::flush;
        
        // Only one every 'outputeverynthtimestep' solutions is output:
        if (timestepindex%outputeverynthtimestep == 0)
        {
            output[0][timestepindex/outputeverynthtimestep] = u;
            output[1][timestepindex/outputeverynthtimestep] = v;
            output[2][timestepindex/outputeverynthtimestep] = a;
        }
        timestepindex++;
    }
    std::cout << std::endl;
    
    // Remove all time derivatives from the universe:
    universe::xdtxdtdtx = {{},{},{}};
    
    return output;
}


