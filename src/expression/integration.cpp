#include "integration.h"


integration::integration(int physreg, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    myexpression = {tointegrate};
    myintegrationorderdelta = integrationorderdelta;
    myblocknumber = blocknumber;
    myphysicalregion = physreg;
}

integration::integration(int physreg, expression meshdeform, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    // Get only the disjoint regions with highest dimension elements:
    std::vector<int> selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions();

    if (not(meshdeform.isharmonicone(selecteddisjregs)))
    {
        logs log;
        log.msg() << "Error in 'integration' object: the mesh deformation expression cannot be multiharmonic (only constant harmonic 1)" << std::endl;
        log.error();
    }

    myexpression = {tointegrate};
    myintegrationorderdelta = integrationorderdelta;
    myblocknumber = blocknumber;
    myphysicalregion = physreg;
    mymeshdeform = {meshdeform};
}

integration::integration(int physreg, int numcoefharms, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    myexpression = {tointegrate};
    myintegrationorderdelta = integrationorderdelta;
    myblocknumber = blocknumber;
    myphysicalregion = physreg;
    mynumcoefharms = numcoefharms;
}

integration::integration(int physreg, int numcoefharms, expression meshdeform, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    // Get only the disjoint regions with highest dimension elements:
    std::vector<int> selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions();

    if (not(meshdeform.isharmonicone(selecteddisjregs)))
    {
        logs log;
        log.msg() << "Error in 'integration' object: the mesh deformation expression cannot be multiharmonic (only constant harmonic 1)" << std::endl;
        log.error();
    }

    myexpression = {tointegrate};
    myintegrationorderdelta = integrationorderdelta;
    myblocknumber = blocknumber;
    myphysicalregion = physreg;
    mymeshdeform = {meshdeform};
    mynumcoefharms = numcoefharms;
}

expression integration::getexpression(void) 
{ 
    return myexpression[0]; 
}

expression integration::getmeshdeform(void) 
{ 
    return mymeshdeform[0]; 
}

void integration::print(void)
{
    std::cout << "Integration of expression" << std::endl;
    myexpression[0].print();
    std::cout << std::endl << "on physical region " << myphysicalregion << std::endl;
    if (mymeshdeform.size() == 0)
        std::cout << "on undeformed mesh" << std::endl;
    else
    {
        std::cout << "on mesh deformed by" << std::endl;
        mymeshdeform[0].print();
        std::cout << std::endl;
    }
    std::cout << "Contribution is added to block " << myblocknumber << std::endl;
    if (mynumcoefharms >= 0)
        std::cout << "An FFT is performed on the coef using " << mynumcoefharms << " time computations" << std::endl;
}



