#include "opdof.h"


void opdof::setspacederivative(int whichderivative)
{
    // Make sure a single space derivative is applied.
    if (spacederivative != 0)
    {
        std::cout << "Error in 'opdof' object: cannot apply more than one space derivative to the dof" << std::endl;
        abort();
    }
    spacederivative = whichderivative;
}

void opdof::increasetimederivativeorder(int amount)
{
    timederivativeorder += amount;

    if (not(myfield->ismultiharmonic()) && timederivativeorder > 2)
    {
        std::cout << "Error in 'opdof' object: time derivative order can exceed 2 only for multiharmonic fields" << std::endl;
        abort();
    }
}

bool opdof::isharmonicone(std::vector<int> disjregs)
{
    std::vector<int> myharms = myfield->getharmonics();
    return (myharms.size() == 1 && myharms[0] == 1);
}

bool opdof::isvalueorientationdependent(std::vector<int> disjregs)
{
    std::cout << "Error in 'opdof' object: 'isvalueorientationdependent' was called (this is not expected)" << std::endl;
    abort();
}

std::shared_ptr<operation> opdof::copy(void)
{
    std::shared_ptr<opdof> op(new opdof(myfield, myphysicalregion));
    *op = *this;
    return op;
}

void opdof::print(void)
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

    std::cout << "dof(";
    myfield->print();
    std::cout << ")";
}
