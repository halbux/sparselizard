#include "coefmanager.h"
#include "universe.h"

coefmanager::coefmanager(std::string fieldtypename)
{
    myfieldtypename = fieldtypename;

    if (myfieldtypename == "x" || myfieldtypename == "y" || myfieldtypename == "z")
    {
        std::cout << "Error in 'coefmanager' object: x, y and z coordinate fields are not supported here" << std::endl;
        abort();
    };
        
    // Preallocate 'coefs' for the number of disjoint regions:
    coefs.resize((universe::mymesh->getdisjointregions())->count());
    // Resize coefs to accomodate an inital order 1 interpolated field:
    for (int i = 0; i < coefs.size(); i++)
        fitinterpolationorder(i, 1);
}

bool coefmanager::isdefined(int disjreg, int formfunctionindex)
{
    return (formfunctionindex < coefs[disjreg].size());
}

void coefmanager::fitinterpolationorder(int disjreg, int interpolationorder)
{
    // Get the element type number in the current disjoint region:
    disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();
    int elementtypenumber = mydisjointregions->getelementtypenumber(disjreg);
    int elementdimension = mydisjointregions->getelementdimension(disjreg);

    // Get the number of form functions associated to dimension elementdimension:
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, myfieldtypename);
    int numberofformfunctions = myformfunction->count(interpolationorder, elementdimension, 0);
    
    if (coefs[disjreg].size() != numberofformfunctions)
        coefs[disjreg].resize(numberofformfunctions);
}

double coefmanager::getcoef(int disjreg, int formfunctionindex, int elementindexindisjointregion)
{
    if (coefs[disjreg][formfunctionindex].size() != 0)
        return coefs[disjreg][formfunctionindex][elementindexindisjointregion];
    else
        return 0;
}

void coefmanager::setcoef(int disjreg, int formfunctionindex, int elementindexindisjointregion, double val)
{
    if (coefs[disjreg][formfunctionindex].size() != 0)
        coefs[disjreg][formfunctionindex][elementindexindisjointregion] = val;
    else
    {
        // This is rarely called and can thus be slower:
        disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();
        int numberofelements = mydisjointregions->countelements(disjreg);
        coefs[disjreg][formfunctionindex].resize(numberofelements); // Filled with zeros.

        coefs[disjreg][formfunctionindex][elementindexindisjointregion] = val;
    }
}
