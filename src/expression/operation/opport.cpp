#include "opport.h"


opport::opport(std::shared_ptr<rawport> portin)
{
    myport = portin;
}

void opport::increasetimederivativeorder(int amount)
{
    timederivativeorder += amount;

    if (timederivativeorder > 2)
    {
        logs log;
        log.msg() << "Error in 'opport' object: time derivative order cannot exceed 2" << std::endl;
        log.error();
    }
}

std::vector<std::vector<densemat>> opport::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    logs log;
    log.msg() << "Error in 'opport' object: cannot interpolate a port" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

densemat opport::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    logs log;
    log.msg() << "Error in 'opport' object: cannot interpolate a port" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

std::shared_ptr<rawport> opport::getportpointer(void)
{
    return myport;
}

bool opport::isharmonicone(std::vector<int> disjregs)
{
    std::vector<int> myharms = myport->getharmonics();
    return (myharms.size() == 1 && myharms[0] == 1);
}

std::shared_ptr<operation> opport::copy(void)
{
    std::shared_ptr<opport> op(new opport(myport));
    *op = *this;
    return op;
}

void opport::print(void)
{
    for (int i = 0; i < timederivativeorder; i++)
        std::cout << "dt";
        
    std::string portname = myport->getname();

    if (portname.size() == 0)
        std::cout << "port";
    else
        std::cout << portname;
}

