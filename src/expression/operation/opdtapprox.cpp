#include "opdtapprox.h"


opdtapprox::opdtapprox(int dtorder, std::shared_ptr<operation> arg, double initdtx, double initdtdtx)
{
    mydtorder = dtorder;
    
    myarg = arg;
    
    mydtx = initdtx;
    mydtdtx = initdtdtx;
}

std::vector<std::vector<densemat>> opdtapprox::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    if (universe::fundamentalfrequency > 0)
    {
        logs log;
        log.msg() << "Error in 'opdtapprox' object: this object cannot be used in the frequency domain" << std::endl;
        log.error();
    }

    densemat output;

    if (mydtorder == 1)
        output = densemat(elemselect.countinselection(), evaluationcoordinates.size()/3, mydtx);
    else
        output = densemat(elemselect.countinselection(), evaluationcoordinates.size()/3, mydtdtx);
        
    // The constant is on the cos0 harmonic:
    return {{},{output}};
}

std::shared_ptr<operation> opdtapprox::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);
    
    return shared_from_this();
}

std::shared_ptr<operation> opdtapprox::copy(void)
{
    std::shared_ptr<opdtapprox> op(new opdtapprox(mydtorder, myarg, mydtx, mydtdtx));
    *op = *this;
    
    universe::opdtapproxes.push_back(op);
            
    return op;
}

void opdtapprox::nextimpliciteuler(double tinit, double dt)
{
    if (mydtorder == 2)
    {
        logs log;
        log.msg() << "Error in 'opdtapprox' object: cannot approximate a second order time derivative with implicit Euler" << std::endl;
        log.error();
    }

    mydtx = (myarg->evaluateattime(tinit+dt) - myarg->evaluateattime(tinit))/dt;
    mydtdtx = 0;
}

void opdtapprox::nextgenalpha(double beta, double gamma, double alphaf, double alpham, double tinit, double dt)
{
    if (mydtbkps.size() == 0)
        mydtbkps = {mydtx, mydtdtx};
 
    double v = mydtbkps[0], a = mydtbkps[1];
    
    double u = myarg->evaluateattime(tinit);
    double unext = myarg->evaluateattime(tinit+dt);

    double anext = 1.0/(beta*dt*dt)*( unext-u - dt*v - dt*dt*(0.5-beta)*a );
    double vnext = v + (dt*(1-gamma))*a + (gamma*dt)*anext;
    
    mydtx = vnext;
    mydtdtx = anext;
}

void opdtapprox::approvetimestep(void)
{
    mydtbkps = {};
}

void opdtapprox::print(void)
{
    if (mydtorder == 1)
        std::cout << "dt(";
    else
        std::cout << "dtdt(";
        
    myarg->print();
    
    std::cout << ")";
}
