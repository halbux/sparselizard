#include "impliciteuler.h"

impliciteuler::impliciteuler(formulation formul, vec dtxinit, int verbosity, std::vector<bool> isrhskcconstant)
{
    myverbosity = verbosity;

    myformulation = formul;
    if (myformulation.ismassmatrixdefined())
    {
        std::cout << "Error in 'impliciteuler' object: formulation provided cannot have a mass matrix (use another time resolution algorithm)" << std::endl;
        abort();  
    }
    
    dtx = dtxinit;
    isconstant = isrhskcconstant;
    
    if (isconstant.size() != 3)
    {
        std::cout << "Error in 'impliciteuler' object: expected a length 3 vector as fourth argument" << std::endl;
        abort();  
    }
}

void impliciteuler::settimederivative(vec sol)
{
    dtx = sol;
}

void impliciteuler::setadaptivity(double tol, double mints, double maxts, double reffact, double coarfact, double coarthres)
{
    if (tol < 0 || mints < 0 || maxts < 0 || reffact < 0 || coarfact < 0 || coarthres < 0)
    {
        std::cout << "Error in 'impliciteuler' object: expected positive arguments for adaptivity" << std::endl;
        abort();  
    }
    if (mints > maxts)
    {
        std::cout << "Error in 'impliciteuler' object: min timestep cannot be larger than max for adaptivity" << std::endl;
        abort();      
    }
    if (reffact > 1)
    {
        std::cout << "Error in 'impliciteuler' object: expected a refinement factor lower than one for adaptivity" << std::endl;
        abort();      
    }
    if (coarfact < 1)
    {
        std::cout << "Error in 'impliciteuler' object: expected a coarsening factor larger than one for adaptivity" << std::endl;
        abort();        
    }
    if (coarthres > 1)
    {
        std::cout << "Error in 'impliciteuler' object: expected a coarsening threshold lower than one for adaptivity" << std::endl;
        abort();        
    }

    mindt = mints; maxdt = maxts; tatol = tol; rfact = reffact; cfact = coarfact; cthres = coarthres;
}

void impliciteuler::presolve(std::vector<formulation> formuls) { tosolvebefore = formuls; }
void impliciteuler::postsolve(std::vector<formulation> formuls) { tosolveafter = formuls; }

void impliciteuler::next(double timestep)
{
    run(true, timestep, -1);
}

int impliciteuler::next(double timestep, int maxnumnlit)
{
    return run(false, timestep, maxnumnlit);
}

int impliciteuler::run(bool islinear, double timestep, int maxnumnlit)
{
    if (timestep < 0 && mindt == -1)
    {
        std::cout << "Error in 'impliciteuler' object: requested an adaptive timestep but adaptivity settings have not been defined" << std::endl;
        abort();
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

    // Get the data from all fields to create the x vector:
    vec x(myformulation);
    x.setdata();
    // Get the initial value of the fields in all other formulations to solve:
    std::vector<vec> presols(tosolvebefore.size()), postsols(tosolveafter.size());
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

    // Time-adaptivity loop:
    int nlit;
    vec xnext, dtxnext;
    while (true)
    {
        // Update and print the time:
        universe::currenttimestep = inittime+dt;
        char spacer = ':';
        if (islinear || myverbosity < 3)
            spacer = ' ';
        if (myverbosity > 1 && istadapt)
            std::cout << "@" << inittime << "+" << dt << "s" << spacer << std::flush;
        if (myverbosity == 1 || myverbosity > 1 && not(istadapt))
            std::cout << "@" << inittime+dt << "s" << spacer << std::flush;
        
        // Make all time derivatives available in the universe:
        universe::xdtxdtdtx = {{},{dtx},{}};
            
        // Nonlinear loop:
        double relchange = 1; nlit = 0;
        xnext = x; dtxnext = dtx;
        while (relchange > nltol && (maxnumnlit <= 0 || nlit < maxnumnlit))
        {
            // Solve all formulations that must be solved at the beginning of the nonlinear loop:
            mathop::solve(tosolvebefore);
            
            vec xtolcalc = xnext;
            
            // Reassemble only the non-constant matrices:
            bool isfirstcall = not(K.isdefined());
            if (isconstant[0] == false || isfirstcall)
            {
                myformulation.generaterhs();
                rhs = myformulation.rhs();
            }
            else
                rhs.updateconstraints();
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
            
            if (islinear == false && myverbosity > 1)
                std::cout << " " << relchange << std::flush;

            nlit++; 
            
            // Solve all formulations that must be solved at the end of the nonlinear loop:
            mathop::solve(tosolveafter);
            
            // Make all time derivatives available in the universe:
            universe::xdtxdtdtx = {{},{dtxnext},{}};
            
            if (islinear)
                break;
        }
        
        if (myverbosity > 2 && islinear == false)
            std::cout << " (" << nlit << "NL it) " << std::flush;
        
        if (istadapt == false)
            break;
        else
        {
            // Deviation from constant dtx to measure the error:
            double errormeasure = dt*(dtxnext - dtx).norm()/xnext.norm();

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
                mathop::setdata(x);
                for (int i = 0; i < presols.size(); i++)
                    mathop::setdata(presols[i]);
                for (int i = 0; i < postsols.size(); i++)
                    mathop::setdata(postsols[i]);
            }
                
            dt = std::min(dt, maxdt);
            dt = std::max(dt, mindt);

            if (istadapt && myverbosity > 2)
                std::cout << "(" << errormeasure << ") " << std::flush;
            
            if (breakit)
                break;
        }
    }
    
    dtx = dtxnext;
    mytimes.push_back(universe::currenttimestep);
    
    return nlit;
}


