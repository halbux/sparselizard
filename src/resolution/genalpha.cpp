#include "genalpha.h"

genalpha::genalpha(formulation formul, vec initdisplacement, vec initspeed, vec initacceleration, std::vector<bool> isrhskcmconstant)
{
    myformulation = formul;
    
    u = initdisplacement;
    v = initspeed;
    a = initacceleration;
    
    if (isrhskcmconstant.size() == 0)
        isconstant = {false,false,false,false};
    else
        isconstant = isrhskcmconstant;
    if (isconstant.size() != 4)
    {
        std::cout << "Error in 'genalpha' object: expected a length 4 or empty vector as fourth argument" << std::endl;
        abort();  
    }
}

void genalpha::setparameter(double rinf)
{
    if (rinf < 0)
    {
        std::cout << "Error in 'genalpha' object: high-frequency dissipation value provided to .setparameter cannot be negative" << std::endl;
        abort();  
    }
    
    // Minimize low frequency dissipation:
    alphaf = rinf/(rinf+1.0);
    alpham = (2.0*rinf-1.0)/(rinf+1.0);
    // Maximize high frequency dissipation:
    beta = 0.25*(1.0-alpham+alphaf)*(1.0-alpham+alphaf);
    // Second order accuracy:
    gamma = 0.5-alpham+alphaf;
}

void genalpha::presolve(std::vector<formulation> formuls) { tosolvebefore = formuls; }
void genalpha::postsolve(std::vector<formulation> formuls) { tosolveafter = formuls; }
        
std::vector<std::vector<vec>> genalpha::runlinear(double starttime, double timestep, double endtime, int outputeverynthtimestep, int verbosity)
{
    return run(true, starttime, timestep, endtime, -1, outputeverynthtimestep, verbosity);
}

std::vector<std::vector<vec>> genalpha::runnonlinear(double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity)
{
    return run(false, starttime, timestep, endtime, maxnumnlit, outputeverynthtimestep, verbosity);
}

std::vector<std::vector<vec>> genalpha::run(bool islinear, double starttime, double timestep, double endtime, int maxnumnlit, int outputeverynthtimestep, int verbosity)
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
    
    vec rhs; mat K, C, M, leftmat, matu, matv, mata;
    
    // Count the number of time steps to step through and the number of vectors to output:
    int numtimesteps = 0; int outputsize = 0;
    for (double t = starttime; t <= endtime; t = t + timestep)
    {
        if (numtimesteps%outputeverynthtimestep == 0)
            outputsize++;
        numtimesteps++;
    }
    
    
    // Start the generalized alpha iteration:
    std::cout << "Generalized alpha (af " << alphaf << ", am " << alpham << ", b " << beta << ", g " << gamma << ") for " << numtimesteps << " timesteps in range " << starttime << " to " << endtime << " sec:" << std::endl;
    std::vector<std::vector<vec>> output(3, std::vector<vec>(outputsize));
    output[0][0] = u; output[1][0] = v; output[2][0] = a;
    
    // We already have everything for time step 0 so we start at 1:
    int timestepindex = 1;
    for (double t = starttime + timestep; t <= endtime; t = t + timestep)
    {        
        std::cout << timestepindex << "@" << t << "sec" << std::flush;

        mathop::settime(t-alphaf*timestep);
        
        // Nonlinear loop:
        double relchange = 1; int nlit = 0;
        vec unext = u, vnext = v, anext = a;
        while (relchange > tol && (maxnumnlit <= 0 || nlit < maxnumnlit))
        {
            // Solve all formulations that must be solved at the beginning of the nonlinear loop:
            mathop::solve(tosolvebefore);


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
                leftmat = (1.0-alpham)*M + (1.0-alphaf) * ( (gamma*timestep)*C + (beta*timestep*timestep)*K );
                leftmat.reuselu();
                
                matu = -K;
                matv = -(1.0-alphaf)*timestep*K-C;
                mata = (1.0-alphaf) * ( ((gamma-1.0)*timestep)*C+((beta-0.5)*timestep*timestep)*K ) - alpham*M;
            }
            
            // Update the acceleration. 
            // The acceleration is imposed on the Dirichlet-constrained dofs.
            // The displacement update relation is used to make sure the 
            // acceleration constraint leads to the exact constrained 
            // displacement at the next time step:
            vec unextdirichlet(myformulation); 
            unextdirichlet.updateconstraints();
            vec anextdirichlet = 1.0/(beta*timestep*timestep)*( unextdirichlet-u - timestep*v - timestep*timestep*(0.5-beta)*a );
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


        	// Solve all formulations that must be solved at the end of the nonlinear loop:
        	mathop::solve(tosolveafter);
        	
            
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


