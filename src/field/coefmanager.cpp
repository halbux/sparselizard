#include "coefmanager.h"

coefmanager::coefmanager(std::string fieldtypename, disjointregions* drs)
{
    myfieldtypename = fieldtypename;
    mydisjointregions = *drs;

    if (myfieldtypename == "x" || myfieldtypename == "y" || myfieldtypename == "z")
    {
        logs log;
        log.msg() << "Error in 'coefmanager' object: x, y and z coordinate fields are not supported here" << std::endl;
        log.error();
    };
        
    // Preallocate 'coefs' for the number of disjoint regions:
    coefs.resize(mydisjointregions.count());
    // Resize coefs to accomodate an inital order 1 interpolated field:
    for (int i = 0; i < coefs.size(); i++)
        fitinterpolationorder(i, 1);
}

bool coefmanager::isdefined(int disjreg, int formfunctionindex)
{
    return (formfunctionindex < coefs[disjreg].size());
}

int coefmanager::countformfunctions(int disjreg)
{
    return coefs[disjreg].size();
}

void coefmanager::fitinterpolationorder(int disjreg, int interpolationorder)
{
    // Get the element type number in the current disjoint region:
    int elementtypenumber = mydisjointregions.getelementtypenumber(disjreg);
    int elementdimension = mydisjointregions.getelementdimension(disjreg);

    // Get the number of form functions associated to dimension elementdimension:
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, myfieldtypename);
    int numberofformfunctions = myformfunction->count(interpolationorder, elementdimension, 0);
    
    if (coefs[disjreg].size() != numberofformfunctions)
        coefs[disjreg].resize(numberofformfunctions);
}

double coefmanager::getcoef(int disjreg, int formfunctionindex, int elementindexindisjointregion)
{
    if (formfunctionindex < coefs[disjreg].size() && coefs[disjreg][formfunctionindex].size() != 0)
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
        int numberofelements = mydisjointregions.countelements(disjreg);
        coefs[disjreg][formfunctionindex].resize(numberofelements); // Filled with zeros.

        coefs[disjreg][formfunctionindex][elementindexindisjointregion] = val;
    }
}

void coefmanager::print(bool databoundsonly)
{
    std::cout << std::endl << "Field of type " << myfieldtypename << ":" << std::endl;
    
    for (int d = 0; d < coefs.size(); d++)
    {
        for (int ff = 0; ff < coefs[d].size(); ff++)
        {
            if (coefs[d][ff].size() == 0)
                continue;
                
            std::cout << std::endl << "--> Disjoint region " << d << ", shape function " << ff << ":" << std::endl;
            
            double datamin = coefs[d][ff][0];
            double datamax = coefs[d][ff][0];
            
            for (int e = 0; e < coefs[d][ff].size(); e++)
            {
                double val = coefs[d][ff][e];
                
                if (databoundsonly == false)
                    std::cout << val << " ";
                if (val < datamin)
                    datamin = val;
                if (val > datamax)
                    datamax = val;
            }
            if (databoundsonly == false)
                std::cout << std::endl;

            std::cout << "Data min/max: " << datamin << " / " << datamax << std::endl;
        }
    }
    std::cout << std::endl;
}
