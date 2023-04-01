#include "optf.h"


void optf::setspacederivative(int whichderivative)
{ 
    // Make sure a single space derivative is applied.
    if (spacederivative != 0)
    {
        logs log;
        log.msg() << "Error in 'optf' object: cannot apply more than one space derivative to the test function" << std::endl;
        log.error();
    }
    spacederivative = whichderivative; 
}

void optf::increasetimederivativeorder(int amount)
{    
    logs log;
    log.msg() << "Error in 'optf' object: cannot apply a time derivative to the test function" << std::endl;
    log.error();
}

bool optf::isharmonicone(std::vector<int> disjregs) 
{ 
    std::vector<int> myharms = myfield->getharmonics(); 
    return (myharms.size() == 1 && myharms[0] == 1);
}

bool optf::isvalueorientationdependent(std::vector<int> disjregs)
{
    logs log;
    log.msg() << "Error in 'optf' object: 'isvalueorientationdependent' was called (this is not expected)" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

std::shared_ptr<operation> optf::copy(void)
{
    std::shared_ptr<optf> op(new optf(myfield, myphysicalregion));
    *op = *this;
    return op;
}

void optf::print(void)
{
    for (int i = 0; i < timederivativeorder; i++)
        std::cout << "dt";
    
    std::vector<std::string> nonedxdydz = {"","dx","dy","dz"};
    std::vector<std::string> nonedkidetadphi = {"","dki","deta","dphi"};
    std::vector<std::string> nonecompxcompycompz = {"compx","compy","compz"};

    std::cout << nonedxdydz[spacederivative];
    std::cout << nonedkidetadphi[kietaphiderivative];
    // For fields without subfields:
    if (fieldcomponent == -1)
    {
        if (myfield->countformfunctioncomponents() > 1)
            std::cout << nonecompxcompycompz[formfunctioncomponent];
    }
    else
        std::cout << nonecompxcompycompz[fieldcomponent];
    
    std::cout << "tf(";
    myfield->print();
    std::cout << ")";
}
