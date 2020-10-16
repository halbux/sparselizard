#include "coordinategroup.h"
#include "universe.h"


coordinategroup::coordinategroup(std::vector<double>& coords)
{
    mynumcoords = coords.size()/3;
    
    noisethreshold = universe::mymesh->getnodes()->getnoisethreshold();
    double meshsize = 0.0;
    for (int i = 0; i < 3; i++)
        meshsize += universe::mymesh->getnodes()->getgeometrydimension(i);
    
    // Define the number of slices in the x, y and z direction:
    double powertouse = universe::mymesh->getmeshdimension();
    int ns = std::pow(mynumcoords, 1.0/powertouse);
    ns = std::max(ns, 1);
    numslices = {ns,ns,ns};
    
    // Get the coordinate x, y and z bounds as well as the distance between slices:
    bounds = myalgorithm::getcoordbounds(coords);
    delta = {std::abs(bounds[1]-bounds[0])/numslices[0], std::abs(bounds[3]-bounds[2])/numslices[1], std::abs(bounds[5]-bounds[4])/numslices[2]};
    
    // This solves the non-existing dimension issues:
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(bounds[2*i+1]-bounds[2*i+0]) < 1e-6*meshsize)
        {
            delta[i] = meshsize;
            numslices[i] = 1;
        }
    }
    
    mygroups = std::shared_ptr<int>(new int[mynumcoords]);
    mygroupcoords = std::shared_ptr<double>(new double[3*mynumcoords]);

    myalgorithm::slicecoordinates(coords, bounds[0], bounds[2], bounds[4], delta[0], delta[1], delta[2], numslices[0], numslices[1], numslices[2], mygroupads, mygroups.get(), mygroupcoords.get());
}

int coordinategroup::countcoordinates(void)
{
    return mynumcoords;
}

void coordinategroup::select(double x, double y, double z, double r)
{
    // Take an extra noise margin to be sure not to miss any candidate slice:
    selx1 = std::floor( ( x-r-noisethreshold[0] - bounds[0] )/delta[0] );
    selx2 = std::floor( ( x+r+noisethreshold[0] - bounds[0] )/delta[0] );
    
    sely1 = std::floor( ( y-r-noisethreshold[1] - bounds[2] )/delta[1] );
    sely2 = std::floor( ( y+r+noisethreshold[1] - bounds[2] )/delta[1] );
    
    selz1 = std::floor( ( z-r-noisethreshold[2] - bounds[4] )/delta[2] );
    selz2 = std::floor( ( z+r+noisethreshold[2] - bounds[4] )/delta[2] );
    
    // Bring in bounds:
    selx1 = std::max(selx1, 0);
    selx1 = std::min(selx1, numslices[0]-1);
    selx2 = std::max(selx2, 0);
    selx2 = std::min(selx2, numslices[0]-1);

    sely1 = std::max(sely1, 0);
    sely1 = std::min(sely1, numslices[1]-1);
    sely2 = std::max(sely2, 0);
    sely2 = std::min(sely2, numslices[1]-1);
    
    selz1 = std::max(selz1, 0);
    selz1 = std::min(selz1, numslices[2]-1);
    selz2 = std::max(selz2, 0);
    selz2 = std::min(selz2, numslices[2]-1);
    
    curselx = selx1; cursely = sely1; curselz = selz1;
    curg = curselx*numslices[1]*numslices[2] + cursely*numslices[2] + curselz;
}

bool coordinategroup::next(void)
{
    if (curselz == selz2)
    {
        if (cursely == sely2)
        {
            if (curselx == selx2)
                return false;
            else
                curselx++;
                
            cursely = sely1;
        }
        else
            cursely++;
        
        curselz = selz1;
    }
    else
        curselz++;
        
    curg = curselx*numslices[1]*numslices[2] + cursely*numslices[2] + curselz;
        
    return true;
}

int coordinategroup::countgroupcoordinates(void)
{
    return (mygroupads[curg+1]-mygroupads[curg]);
}
        
int* coordinategroup::getgroupindexes(void)
{
    return mygroups.get()+mygroupads[curg];
}

double* coordinategroup::getgroupcoordinates(void)
{
    return mygroupcoords.get()+3*mygroupads[curg];
}

