#include "genalpha.h"

genalpha::genalpha(formulation formul, vec initspeed, vec initacceleration, std::vector<bool> isrhskcmconstant)
{
    myformulation = formul;
    
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
    
    alphaf = rinf/(rinf+1.0);
    alpham = (2.0*rinf-1.0)/(rinf+1.0);
    beta = 0.25*(1.0-alpham+alphaf)*(1.0-alpham+alphaf);
    gamma = 0.5-alpham+alphaf;
}

void genalpha::settimederivative(std::vector<vec> sol)
{
    if (sol.size() != 2)
    {
        std::cout << "Error in 'genalpha' object: expected a vector of length two to set the time derivatives" << std::endl;
        abort();  
    }
    v = sol[0]; a = sol[1];
}

void genalpha::presolve(std::vector<formulation> formuls) { tosolvebefore = formuls; }
void genalpha::postsolve(std::vector<formulation> formuls) { tosolveafter = formuls; }
        
void genalpha::runlinear(double timestep, int verbosity, bool autoadvancetime)
{
    run(true, timestep, -1, verbosity, autoadvancetime);
}

int genalpha::runnonlinear(double timestep, int maxnumnlit, int verbosity, bool autoadvancetime)
{
    return run(false, timestep, maxnumnlit, verbosity, autoadvancetime);
}

int genalpha::run(bool islinear, double timestep, int maxnumnlit, int verbosity, bool autoadvancetime)
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
    universe::xdtxdtdtx = {{},{v},{a}};
        
    // Get the data from all fields to create the u vector:
    vec u(myformulation);
    u.setdata();
    
    // Nonlinear loop:
    double relchange = 1; int nlit = 0;
    vec unext = u, vnext = v, anext = a;
    while (relchange > nltol && (maxnumnlit <= 0 || nlit < maxnumnlit))
    {
        // Solve all formulations that must be solved at the beginning of the nonlinear loop:
        mathop::solve(tosolvebefore);

        // Make all time derivatives available in the universe:
        universe::xdtxdtdtx = {{},{vnext},{anext}};
    
        vec utolcalc = unext;
        
        // Reassemble only the non-constant matrices:
        bool isfirstcall = not(K.isdefined());
        if (isconstant[1] == false || isfirstcall)
        {
            myformulation.generatestiffnessmatrix();
            K = myformulation.K(false, true);
        }
        if (isconstant[2] == false || isfirstcall)
        {
            myformulation.generatedampingmatrix();
            C = myformulation.C(false, true);
        }
        if (isconstant[3] == false || isfirstcall)
        {
            myformulation.generatemassmatrix();
            M = myformulation.M(false, false);
        }
        if (isconstant[0] == false || isfirstcall)
        {
            myformulation.generaterhs();
            rhs = myformulation.rhs();
        }
        else
            rhs.updateconstraints();
        
        // Reuse matrices when possible (including the factorization):
        if (isconstant[1] == false || isconstant[2] == false || isconstant[3] == false || isfirstcall || defdt != dt || defbeta != beta || defgamma != gamma || defalphaf != alphaf || defalpham != alpham)
        {
            leftmat = (1.0-alpham)*M + ((1.0-alphaf)*gamma*dt)*C + ((1.0-alphaf)*beta*dt*dt)*K;
            leftmat.reusefactorization();
            
            matu = -K;
            matv = ((alphaf-1.0)*dt)*K-C;
            mata = ((1.0-alphaf)*(gamma-1.0)*dt)*C+((1.0-alphaf)*(beta-0.5)*dt*dt)*K - alpham*M;
            
            defdt = dt; defbeta = beta; defgamma = gamma; defalphaf = alphaf; defalpham = alpham;
        }
        
        // Update the acceleration. 
        // The acceleration is imposed on the Dirichlet-constrained dofs.
        // The displacement update relation is used to make sure the 
        // acceleration constraint leads to the exact constrained 
        // displacement at the next time step:
        vec unextdirichlet(myformulation); 
        unextdirichlet.updateconstraints();
        vec anextdirichlet = (1.0-alpham)/(beta*dt*dt)*( unextdirichlet-u - dt*v - dt*dt*(0.5-beta)*a );
        // Here are the constrained values of the next acceleration:
        intdensematrix constraintindexes = myformulation.getdofmanager()->getconstrainedindexes();
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
        unext = u + dt*v + ((0.5-beta)*dt*dt)*a + (beta*dt*dt)*anext;
        vnext = v + (dt*(1-gamma))*a + (gamma*dt)*anext;
        
        // Update all fields in the formulation:
        mathop::setdata(unext);
        
        relchange = (unext-utolcalc).norm()/unext.norm();
        
        if (islinear == false && verbosity > 1)
            std::cout << " " << relchange << std::flush;

        nlit++; 

        // Solve all formulations that must be solved at the end of the nonlinear loop:
        mathop::solve(tosolveafter);
        
        if (islinear)
            break;
    }

    v = vnext; a = anext;
    
    if (verbosity > 1 && islinear == false)
        std::cout << " (" << nlit << "NL it) " << std::flush;
    
    if (autoadvancetime == false)
        universe::currenttimestep = inittime;
        
    return nlit;
}


