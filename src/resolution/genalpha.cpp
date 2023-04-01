#include "genalpha.h"

genalpha::genalpha(formulation formul, vec initspeed, vec initacceleration, int verbosity, std::vector<bool> isrhskcmconstant)
{
    myverbosity = verbosity;

    myformulation = formul;
    
    v = initspeed;
    a = initacceleration;
    isconstant = isrhskcmconstant;
        
    if (isconstant.size() != 4)
    {
        logs log;
        log.msg() << "Error in 'genalpha' object: expected a length 4 vector as fifth argument" << std::endl;
        log.error();  
    }
}

void genalpha::setparameter(double rinf)
{
    if (rinf < -1e-8)
    {
        logs log;
        log.msg() << "Error in 'genalpha' object: high-frequency dissipation value provided to .setparameter cannot be negative" << std::endl;
        log.error();  
    }
    if (rinf > 1+1e-8)
    {
        logs log;
        log.msg() << "Error in 'genalpha' object: high-frequency dissipation value provided to .setparameter cannot be larger than one" << std::endl;
        log.error();  
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
        logs log;
        log.msg() << "Error in 'genalpha' object: expected a vector of length two to set the time derivatives" << std::endl;
        log.error();  
    }
    v = sol[0]; a = sol[1];
}

void genalpha::setadaptivity(double tol, double mints, double maxts, double reffact, double coarfact, double coarthres)
{
    if (tol < 0 || mints < 0 || maxts < 0 || reffact < 0 || coarfact < 0 || coarthres < 0)
    {
        logs log;
        log.msg() << "Error in 'genalpha' object: expected positive arguments for adaptivity" << std::endl;
        log.error();  
    }
    if (mints > maxts)
    {
        logs log;
        log.msg() << "Error in 'genalpha' object: min timestep cannot be larger than max for adaptivity" << std::endl;
        log.error();      
    }
    if (reffact > 1)
    {
        logs log;
        log.msg() << "Error in 'genalpha' object: expected a refinement factor lower than one for adaptivity" << std::endl;
        log.error();      
    }
    if (coarfact < 1)
    {
        logs log;
        log.msg() << "Error in 'genalpha' object: expected a coarsening factor larger than one for adaptivity" << std::endl;
        log.error();        
    }
    if (coarthres > 1)
    {
        logs log;
        log.msg() << "Error in 'genalpha' object: expected a coarsening threshold lower than one for adaptivity" << std::endl;
        log.error();        
    }

    mindt = mints; maxdt = maxts; tatol = tol; rfact = reffact; cfact = coarfact; cthres = coarthres;
}

void genalpha::presolve(std::vector<formulation> formuls) { tosolvebefore = formuls; }
void genalpha::postsolve(std::vector<formulation> formuls) { tosolveafter = formuls; }
        
void genalpha::next(double timestep)
{
    run(true, timestep, -1);
}

int genalpha::next(double timestep, int maxnumnlit)
{
    return run(false, timestep, maxnumnlit);
}

int genalpha::run(bool islinear, double timestep, int maxnumnlit)
{
    if (timestep < 0 && mindt == -1)
    {
        logs log;
        log.msg() << "Error in 'genalpha' object: requested an adaptive timestep but adaptivity settings have not been defined" << std::endl;
        log.error();
    }

    double inittime = universe::currenttimestep;

    // Adaptive timestep:
    bool istadapt = false;
    if (timestep < 0)
    {
        istadapt = true;
        if (dt < 0)
            dt = mindt;
    }
    else
        dt = timestep;

    // Get the data from all fields to create the u vector:
    vec u(myformulation);
    u.setdata();
    // Get the initial value of the fields in all other formulations to solve:
    std::vector<vec> presols(tosolvebefore.size()), postsols(tosolveafter.size());
    if (istadapt)
    {
        for (int i = 0; i < presols.size(); i++)
        {
            presols[i] = vec(tosolvebefore[i]);
            presols[i].setdata();
        }
        for (int i = 0; i < postsols.size(); i++)
        {
            postsols[i] = vec(tosolveafter[i]);
            postsols[i].setdata();
        }
    }

    // Time-adaptivity loop:
    int nlit;
    vec unext, vnext, anext;
    while (true)
    {
        // Update and print the time:
        universe::currenttimestep = inittime+dt;
        
        // Update all opdtapproxes:
        std::vector<std::shared_ptr<operation>> dtapproxes = universe::getdtapproxes();
        for (int i = 0; i < dtapproxes.size(); i++)
            dtapproxes[i]->nextgenalpha(beta, gamma, alphaf, alpham, inittime, dt);

        if (myverbosity > 1 && istadapt)
            std::cout << "@" << inittime << "+" << dt << "s " << std::flush;
        if (myverbosity > 1 && not(istadapt))
            std::cout << "@" << inittime+dt << "s " << std::flush;
    
        // Make all time derivatives available in the universe:
        universe::xdtxdtdtx = {{},{v},{a}};
        
        // Nonlinear loop:
        double relchange = 1; nlit = 0;
        unext = u; vnext = v; anext = a;
        while (relchange > nltol && (maxnumnlit <= 0 || nlit < maxnumnlit))
        {
            // Solve all formulations that must be solved at the beginning of the nonlinear loop:
            for (int i = 0; i < tosolvebefore.size(); i++)
                tosolvebefore[i].solve();
        
            vec utolcalc = unext;
            
            // Reassemble only the non-constant matrices:
            bool isfirstcall = not(K.isdefined());
            
            universe::currenttimestep = inittime+(1.0-alphaf)*dt;
            if (isconstant[0] == false || isfirstcall)
            {
                myformulation.generaterhs();
                rhs = myformulation.rhs();
            }
            else
                rhs.updateconstraints();
            universe::currenttimestep = inittime+dt;
            
            if (isconstant[1] == false || isfirstcall)
            {
                myformulation.generatestiffnessmatrix();
                K = myformulation.K(false);
            }
            if (isconstant[2] == false || isfirstcall)
            {
                myformulation.generatedampingmatrix();
                C = myformulation.C(false);
            }
            if (isconstant[3] == false || isfirstcall)
            {
                myformulation.generatemassmatrix();
                M = myformulation.M(false);
            }
            
            // Reuse matrices when possible (including the factorization):
            if (isconstant[1] == false || isconstant[2] == false || isconstant[3] == false || isfirstcall || defdt != dt || defbeta != beta || defgamma != gamma || defalphaf != alphaf || defalpham != alpham)
            {
                leftmat = ((1.0-alpham)/(beta*dt*dt))*M + ((1.0-alphaf)*gamma/(beta*dt))*C + (1.0-alphaf)*K;
                leftmat.reusefactorization();
                
                defdt = dt; defbeta = beta; defgamma = gamma; defalphaf = alphaf; defalpham = alpham;
            }
            
            vec vm = ((alpham-1.0)/(beta*dt*dt))*u + ((alpham-1.0)/(beta*dt))*v + ((alpham-1.0)*(0.5-beta)/beta + alpham)*a;
            vec vc = ((alphaf-1.0)*gamma/(beta*dt))*u + (1.0+(alphaf-1.0)*gamma/beta)*v + ((1.0-alphaf)*dt*(1.0-gamma)+(alphaf-1.0)*gamma*dt*(0.5-beta)/beta)*a;
            vec vk = alphaf*u;
            
            // Here are the constrained values of the next solution:
            indexmat constraintindexes = myformulation.getdofmanager()->getconstrainedindexes();
            densemat unextdirichletval = rhs.getpointer()->getvalues(constraintindexes);
            vec rightvec = rhs - (M*vm + C*vc + K*vk);
            // Force the solution on the constrained dofs:
            rightvec.getpointer()->setvalues(constraintindexes, unextdirichletval);
            
            // Update the solution unext.
            unext = relaxationfactor * sl::solve(leftmat, rightvec) + (1.0-relaxationfactor)*unext;
            
            // Update anext and vnext:
            anext = 1.0/(beta*dt*dt)*( unext-u - dt*v - dt*dt*(0.5-beta)*a );
            vnext = v + (dt*(1-gamma))*a + (gamma*dt)*anext;
            
            // Update all fields in the formulation:
            sl::setdata(unext);
            
            relchange = (unext-utolcalc).norm()/unext.norm();
            
            if (islinear == false && myverbosity > 2)
                std::cout << relchange << " " << std::flush;

            nlit++; 

            // Solve all formulations that must be solved at the end of the nonlinear loop:
            for (int i = 0; i < tosolveafter.size(); i++)
                tosolveafter[i].solve();
            
            // Make all time derivatives available in the universe:
            universe::xdtxdtdtx = {{},{vnext},{anext}};
            
            if (islinear)
                break;
        }
        
        if (myverbosity > 1 && islinear == false)
            std::cout << "(" << nlit << "NL it) " << std::flush;
        
        if (istadapt == false)
            break;
        else
        {
            // Deviation from constant v to measure the error:
            double errormeasure = dt*(vnext - v).norm()/unext.norm();

            bool breakit = false;
            if (dt <= mindt || errormeasure <= tatol && (islinear || maxnumnlit <= 0 || nlit < maxnumnlit))
            {
                // If the error is low enough to coarsen the timestep:
                if (errormeasure <= cthres*tatol && (islinear || maxnumnlit <= 0 || nlit < maxnumnlit))
                    dt *= cfact;
                breakit = true;
            }
            else
            {
                dt *= rfact;
                // Reset fields for a new resolution:
                sl::setdata(u);
                for (int i = 0; i < presols.size(); i++)
                    sl::setdata(presols[i]);
                for (int i = 0; i < postsols.size(); i++)
                    sl::setdata(postsols[i]);
            }
                
            dt = std::min(dt, maxdt);
            dt = std::max(dt, mindt);

            if (myverbosity > 2)
                std::cout << "(" << errormeasure << ") " << std::flush;
            
            if (breakit)
                break;
        }
    }
    
    // Approve all opdtapproxes:
    std::vector<std::shared_ptr<operation>> dtapproxes = universe::getdtapproxes();
    for (int i = 0; i < dtapproxes.size(); i++)
        dtapproxes[i]->approvetimestep();
    
    if (myverbosity == 1)
        std::cout << "@" << universe::currenttimestep << "s " << std::flush;
    
    v = vnext; a = anext;
    mytimes.push_back(universe::currenttimestep);
    
    return nlit;
}

