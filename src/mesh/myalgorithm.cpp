#include "myalgorithm.h"

#if defined(__linux__)
#include <parallel/algorithm>
#endif

void myalgorithm::stablecoordinatesort(std::vector<double> noisethreshold, std::vector<double>& coordinates, std::vector<int>& reorderingvector)
{
    // There is a x, y and z coordinate for every node:
    int numberofnodes = coordinates.size()/3;
    
    // 'reorderingvector' gives the relation between the indexes before and after node sorting:
    if (reorderingvector.size() != numberofnodes)
        reorderingvector.resize(numberofnodes);
    // Set 'reorderingvector' to [0 1 2 ...]:
    std::iota(reorderingvector.begin(), reorderingvector.end(), 0);
    // Sort 'reorderingvector' according to 'coordinates' with x > y > z priority order:
    // The < operator is overloaded by a lambda function.
    #if defined(__linux__)
    __gnu_parallel::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
    #else
    std::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
    #endif
        { 
            // First sort according to the x coordinate:
            if (coordinates[elem1*3+0] < coordinates[elem2*3+0] - noisethreshold[0])
                return true;
            if (coordinates[elem1*3+0] > coordinates[elem2*3+0] + noisethreshold[0])
                return false;
            // If it cannot be sorted according to the x coordinate sort according to y:
            if (coordinates[elem1*3+1] < coordinates[elem2*3+1] - noisethreshold[1])
                return true;
            if (coordinates[elem1*3+1] > coordinates[elem2*3+1] + noisethreshold[1])
                return false;
            // Otherwise sort according to z:
            if (coordinates[elem1*3+2] < coordinates[elem2*3+2] - noisethreshold[2])
                return true;
            if (coordinates[elem1*3+2] > coordinates[elem2*3+2] + noisethreshold[2])
                return false;
            // For identical entries make a COHERENT decision for a stable sorting.
            return (elem1 < elem2);
        });
}

void myalgorithm::stablecoordinatesort(std::vector<double> noisethreshold, std::vector<int>& elems, std::vector<double>& coordinates, std::vector<int>& reorderingvector)
{
    int numberofnodes = elems.size();
    
    // 'reorderingvector' gives the relation between the indexes before and after node sorting:
    if (reorderingvector.size() != numberofnodes)
        reorderingvector.resize(numberofnodes);
    // Set 'reorderingvector' to [0 1 2 ...]:
    std::iota(reorderingvector.begin(), reorderingvector.end(), 0);
    // Sort 'reorderingvector' according to 'coordinates' with x > y > z priority order:
    // The < operator is overloaded by a lambda function.
    std::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
        { 
            // First sort according to the integer vector:
            if (elems[elem1] < elems[elem2])
                return true;
            if (elems[elem1] > elems[elem2])
                return false;
            
            // Otherwise sort according to the x coordinate:
            if (coordinates[elem1*3+0] < coordinates[elem2*3+0] - noisethreshold[0])
                return true;
            if (coordinates[elem1*3+0] > coordinates[elem2*3+0] + noisethreshold[0])
                return false;
            // If it cannot be sorted according to the x coordinate sort according to y:
            if (coordinates[elem1*3+1] < coordinates[elem2*3+1] - noisethreshold[1])
                return true;
            if (coordinates[elem1*3+1] > coordinates[elem2*3+1] + noisethreshold[1])
                return false;
            // Otherwise sort according to z:
            if (coordinates[elem1*3+2] < coordinates[elem2*3+2] - noisethreshold[2])
                return true;
            if (coordinates[elem1*3+2] > coordinates[elem2*3+2] + noisethreshold[2])
                return false;
            // For identical entries make a COHERENT decision for a stable sorting.
            return (elem1 < elem2);
        });
}

int myalgorithm::removeduplicatedcoordinates(std::vector<double> noisethreshold, std::vector<double>& coordinates, std::vector<int>& renumberingvector)
{
    // There is a x, y and z coordinate for every nodes:
    int numberofnodes = coordinates.size()/3;
    
    if (renumberingvector.size() != numberofnodes)
        renumberingvector.resize(numberofnodes);
    
    if (numberofnodes == 0)
        return 0;
    
    int newnodenumber = 0;
    renumberingvector[0] = 0;
    for (int i = 1; i < numberofnodes; i++)
    {
        // If the node is not close enough to the previous one (i.e. is not a duplicate):
        if (std::abs(coordinates[3*i+0] - coordinates[3*(i-1)+0]) > noisethreshold[0] || std::abs(coordinates[3*i+1] - coordinates[3*(i-1)+1]) > noisethreshold[1] || std::abs(coordinates[3*i+2] - coordinates[3*(i-1)+2]) > noisethreshold[2])
            newnodenumber++;

        // Create a vector whose ith node gives the new node number for node i:
        renumberingvector[i] = newnodenumber;
    }

    return newnodenumber + 1;
}

void myalgorithm::stablesort(std::vector<int>& tosort, std::vector<int>& reorderingvector)
{
    if (reorderingvector.size() != tosort.size())
        reorderingvector.resize(tosort.size());
    
    // Set 'reorderingvector' to [0 1 2 ...]:
    std::iota(reorderingvector.begin(), reorderingvector.end(), 0);
    // Sort 'reorderingvector' according to 'tosort':
    // The < operator is overloaded by a lambda function.
    std::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
        { 
            if (tosort[elem1] < tosort[elem2])
                return true;
            if (tosort[elem1] > tosort[elem2])
                return false;
            // For identical entries make a COHERENT decision for a stable sorting.
            return elem1 < elem2;
        });
}

void myalgorithm::stablesort(double noisethreshold, std::vector<double>& tosort, std::vector<int>& reorderingvector)
{
    if (reorderingvector.size() != tosort.size())
        reorderingvector.resize(tosort.size());
    
    // Set 'reorderingvector' to [0 1 2 ...]:
    std::iota(reorderingvector.begin(), reorderingvector.end(), 0);
    // Sort 'reorderingvector' according to 'tosort':
    // The < operator is overloaded by a lambda function.
    std::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
        { 
            if (tosort[elem1] < tosort[elem2] - noisethreshold)
                return true;
            if (tosort[elem1] > tosort[elem2] + noisethreshold)
                return false;
            // For identical entries make a COHERENT decision for a stable sorting.
            return elem1 < elem2;
        });
}

void myalgorithm::stablesort(double noisethreshold, std::vector<double>& tosort, std::vector<int>& reorderingvector, int blocklen)
{
    if (reorderingvector.size() != tosort.size()/blocklen)
        reorderingvector.resize(tosort.size()/blocklen);
    
    // Set 'reorderingvector' to [0 1 2 ...]:
    std::iota(reorderingvector.begin(), reorderingvector.end(), 0);
    // Sort 'reorderingvector' according to 'tosort':
    // The < operator is overloaded by a lambda function.
    std::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
        { 
            for (int i = 0; i < blocklen; i++)
            {
                if (tosort[elem1*blocklen+i] < tosort[elem2*blocklen+i] - noisethreshold)
                    return true;
                if (tosort[elem1*blocklen+i] > tosort[elem2*blocklen+i] + noisethreshold)
                    return false;
            }
            // For identical entries make a COHERENT decision for a stable sorting.
            return elem1 < elem2;
        });
}


bool sortfun(const std::tuple<int,int,double>& elem1, const std::tuple<int,int,double>& elem2)
{ 
    if (std::get<0>(elem1) < std::get<0>(elem2))
        return true;
    if (std::get<0>(elem1) > std::get<0>(elem2))
        return false;
    if (std::get<1>(elem1) < std::get<1>(elem2))
        return true;
    if (std::get<1>(elem1) >= std::get<1>(elem2))
        return false;
}
    
void myalgorithm::tuple3sort(std::vector<std::tuple<int,int,double>>& tosort)
{
    // Parallel sort on Linux only for now:
    #if defined(__linux__)
    __gnu_parallel::sort(tosort.begin(), tosort.end(), sortfun);
    #else
    std::sort(tosort.begin(), tosort.end(), sortfun);
    #endif
}

void myalgorithm::slicecoordinates(double noisethreshold, std::vector<double>& toslice, double minval, double delta, int numslices, std::vector<std::vector<int>>& slices)
{
    int num = toslice.size();
    
    // To be sure not to miss any value:
    minval -= noisethreshold;
    delta += 2.0*noisethreshold;
    
    slices = {};
    
    // Create a vector giving the slice in which each value is:
    std::vector<int> inslice(num);
    std::vector<int> countinslice(numslices,0);
    for (int i = 0; i < num; i++)
    {
        int curslice = std::floor((toslice[i]-minval)/delta);
        if (curslice >= 0 && curslice < numslices)
        {
            inslice[i] = curslice;
            countinslice[curslice]++;
        }
        else
        {
            std::cout << "Error in 'myalgorithm': value to slice is out of range (numerical noise threshold exceeded)" << std::endl;
            abort();
        }
    }
    
    // Preallocate 'slices':
    for (int s = 0; s < numslices; s++)
        slices.push_back(std::vector<int>(countinslice[s]));
    // Populate 'slices':
    std::vector<int> curposinslice(numslices,0);
    for (int i = 0; i < num; i++)
    {
        int s = inslice[i];
        slices[s][curposinslice[s]] = i;
        curposinslice[s]++;
    }
}

std::vector<double> myalgorithm::getcoordbounds(std::vector<double>& coords)
{
    int numcoords = coords.size()/3;
    
    if (numcoords == 0)
        return {};
    
    std::vector<double> output = {coords[0],coords[0], coords[1],coords[1], coords[2],coords[2]};
    
    for (int i = 1; i < numcoords; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            // Min:
            if (output[2*j+0] <= coords[3*i+j]) {}
            else
                output[2*j+0] = coords[3*i+j];
            // Max:
            if (output[2*j+1] >= coords[3*i+j]) {}
            else
                output[2*j+1] = coords[3*i+j];
        }
    }
    
    return output;
}

std::vector<int> myalgorithm::intersect(std::vector<int> a, std::vector<int> b)
{
    // Sort both vectors:
    std::sort(a.begin(), a.end()); 
    std::sort(b.begin(), b.end());
    
    std::vector<int> output(a.size()+b.size());
    std::vector<int>::iterator it;
    
    // Intersect:
    it = std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), output.begin());
    output.resize(it-output.begin()); 
    
    return output;
}

void myalgorithm::csrtoijk(int numberofrows, int* csrrows, int* ijkrows)
{
    // Loop on all rows:
    int index = 0;
    for (int i = 0; i < numberofrows; i++)
    {
        // Loop on all columns:
        for (int j = csrrows[i]; j < csrrows[i+1]; j++)
        {
            ijkrows[index] = i;
            index++;
        }
    }
}

int myalgorithm::getroot(std::vector<polynomial>& poly, std::vector<double>& rhs, std::vector<double>& guess, double boxsize, double tol, int maxit)
{
    int it = 0;
    double deltaki = 1.0, deltaeta = 1.0, deltaphi = 1.0;
    
    // Limit the jumps to a fraction of the box size:
    double maxjump = 0.45;

    switch (poly.size())
    {
        case 1:
        {
            while (std::abs(deltaki) > tol)
            {
                if (it < maxit)
                {
                    double f = poly[0].evalat(guess, 0)[0] - rhs[0];
                    double jac11 = poly[0].evalat(guess, 1)[0];
                    
                    deltaki = -1.0/jac11 * f;
                    
                    double scaling = maxjump*boxsize / ( std::abs(deltaki) );
                    if (scaling < 1.0)
                        deltaki *= scaling;
                    
                    guess[0] += deltaki;
                    
                    if (std::abs(guess[0]) > boxsize)
                        return 0;
                }
                else
                    return -1;
                it++;
            }
            break;
        }
        case 2:
        {
            while (std::abs(deltaki) > tol || std::abs(deltaeta) > tol)
            {
                if (it < maxit)
                {
                    double f = poly[0].evalat(guess, 0)[0] - rhs[0];
                    double g = poly[1].evalat(guess, 0)[0] - rhs[1];
                    
                    double jac11 = poly[0].evalat(guess, 1)[0];
                    double jac12 = poly[0].evalat(guess, 2)[0];
                    double jac21 = poly[1].evalat(guess, 1)[0];
                    double jac22 = poly[1].evalat(guess, 2)[0];
                    
                    double invdet = 1.0/(jac11*jac22 - jac12*jac21);
                    
                    deltaki = -invdet * (jac22*f - jac12*g);
                    deltaeta = -invdet * (-jac21*f + jac11*g);
                    
                    double scaling = maxjump*boxsize / ( std::abs(deltaki)+std::abs(deltaeta) );
                    if (scaling < 1.0)
                    {
                        deltaki *= scaling;
                        deltaeta *= scaling;
                    }
                    
                    guess[0] += deltaki;
                    guess[1] += deltaeta;
                    
                    if (std::abs(guess[0]) > boxsize || std::abs(guess[1]) > boxsize)
                        return 0;
                }
                else
                    return -1;
                it++;
            }
            break;
        }
        case 3:
        {
            while (std::abs(deltaki) > tol || std::abs(deltaeta) > tol || std::abs(deltaphi) > tol)
            {
                if (it < maxit)
                {
                    double f = poly[0].evalat(guess, 0)[0] - rhs[0];
                    double g = poly[1].evalat(guess, 0)[0] - rhs[1];
                    double h = poly[2].evalat(guess, 0)[0] - rhs[2];
                    
                    double jac11 = poly[0].evalat(guess, 1)[0];
                    double jac12 = poly[0].evalat(guess, 2)[0];
                    double jac13 = poly[0].evalat(guess, 3)[0];
                    double jac21 = poly[1].evalat(guess, 1)[0];
                    double jac22 = poly[1].evalat(guess, 2)[0];
                    double jac23 = poly[1].evalat(guess, 3)[0];
                    double jac31 = poly[2].evalat(guess, 1)[0];
                    double jac32 = poly[2].evalat(guess, 2)[0];
                    double jac33 = poly[2].evalat(guess, 3)[0];
                    
                    double invdet = 1.0/(jac11*jac22*jac33 - jac11*jac23*jac32 - jac12*jac21*jac33 + jac12*jac23*jac31 + jac13*jac21*jac32 - jac13*jac22*jac31);
                    
                    deltaki = -invdet * ((jac22*jac33 - jac23*jac32)*f - (jac12*jac33 - jac13*jac32)*g + (jac12*jac23 - jac13*jac22)*h);
                    deltaeta = -invdet * (-(jac21*jac33 - jac23*jac31)*f + (jac11*jac33 - jac13*jac31)*g - (jac11*jac23 - jac13*jac21)*h);
                    deltaphi = -invdet * ((jac21*jac32 - jac22*jac31)*f - (jac11*jac32 - jac12*jac31)*g + (jac11*jac22 - jac12*jac21)*h);
                    
                    double scaling = maxjump*boxsize / ( std::abs(deltaki)+std::abs(deltaeta)+std::abs(deltaphi) );
                    if (scaling < 1.0)
                    {
                        deltaki *= scaling;
                        deltaeta *= scaling;
                        deltaphi *= scaling;
                    }
                    
                    guess[0] += deltaki;
                    guess[1] += deltaeta;
                    guess[2] += deltaphi;
                    
                    if (std::abs(guess[0]) > boxsize || std::abs(guess[1]) > boxsize || std::abs(guess[2]) > boxsize)
                        return 0;
                }
                else
                    return -1;
                it++;
            }
            break;
        }
    }
    return 1;
}

int myalgorithm::getroot(polynomials& polys, std::vector<double>& rhs, std::vector<double>& guess, double boxsize, double tol, int maxit)
{
    int it = 0;
    double deltaki = 1.0, deltaeta = 1.0, deltaphi = 1.0;
    
    // Limit the jumps to a fraction of the box size:
    double maxjump = 0.45;

    switch (polys.count())
    {
        case 2:
        {
            while (std::abs(deltaki) > tol)
            {
                if (it < maxit)
                {
                    std::vector<double> evaled;
                    polys.evalatsingle(guess, evaled);
                
                    double f = evaled[0] - rhs[0];
                    double jac11 = evaled[1];
                    
                    deltaki = -1.0/jac11 * f;
                    
                    double scaling = maxjump*boxsize / ( std::abs(deltaki) );
                    if (scaling < 1.0)
                        deltaki *= scaling;
                    
                    guess[0] += deltaki;
                    
                    if (std::abs(guess[0]) > boxsize)
                        return 0;
                }
                else
                    return -1;
                it++;
            }
            break;
        }
        case 6:
        {
            while (std::abs(deltaki) > tol || std::abs(deltaeta) > tol)
            {
                if (it < maxit)
                {
                    std::vector<double> evaled;
                    polys.evalatsingle(guess, evaled);
                    
                    double f = evaled[0] - rhs[0];
                    double g = evaled[3] - rhs[1];
                    
                    double jac11 = evaled[1];
                    double jac12 = evaled[2];
                    double jac21 = evaled[4];
                    double jac22 = evaled[5];
                    
                    double invdet = 1.0/(jac11*jac22 - jac12*jac21);
                    
                    deltaki = -invdet * (jac22*f - jac12*g);
                    deltaeta = -invdet * (-jac21*f + jac11*g);
                    
                    double scaling = maxjump*boxsize / ( std::abs(deltaki)+std::abs(deltaeta) );
                    if (scaling < 1.0)
                    {
                        deltaki *= scaling;
                        deltaeta *= scaling;
                    }
                    
                    guess[0] += deltaki;
                    guess[1] += deltaeta;
                    
                    if (std::abs(guess[0]) > boxsize || std::abs(guess[1]) > boxsize)
                        return 0;
                }
                else
                    return -1;
                it++;
            }
            break;
        }
        case 12:
        {
            while (std::abs(deltaki) > tol || std::abs(deltaeta) > tol || std::abs(deltaphi) > tol)
            {
                if (it < maxit)
                {
                    std::vector<double> evaled;
                    polys.evalatsingle(guess, evaled);
                    
                    double f = evaled[0] - rhs[0];
                    double g = evaled[4] - rhs[1];
                    double h = evaled[8] - rhs[2];
                    
                    double jac11 = evaled[1];
                    double jac12 = evaled[2];
                    double jac13 = evaled[3];
                    double jac21 = evaled[5];
                    double jac22 = evaled[6];
                    double jac23 = evaled[7];
                    double jac31 = evaled[9];
                    double jac32 = evaled[10];
                    double jac33 = evaled[11];
                    
                    double invdet = 1.0/(jac11*jac22*jac33 - jac11*jac23*jac32 - jac12*jac21*jac33 + jac12*jac23*jac31 + jac13*jac21*jac32 - jac13*jac22*jac31);
                    
                    deltaki = -invdet * ((jac22*jac33 - jac23*jac32)*f - (jac12*jac33 - jac13*jac32)*g + (jac12*jac23 - jac13*jac22)*h);
                    deltaeta = -invdet * (-(jac21*jac33 - jac23*jac31)*f + (jac11*jac33 - jac13*jac31)*g - (jac11*jac23 - jac13*jac21)*h);
                    deltaphi = -invdet * ((jac21*jac32 - jac22*jac31)*f - (jac11*jac32 - jac12*jac31)*g + (jac11*jac22 - jac12*jac21)*h);
                    
                    double scaling = maxjump*boxsize / ( std::abs(deltaki)+std::abs(deltaeta)+std::abs(deltaphi) );
                    if (scaling < 1.0)
                    {
                        deltaki *= scaling;
                        deltaeta *= scaling;
                        deltaphi *= scaling;
                    }
                    
                    guess[0] += deltaki;
                    guess[1] += deltaeta;
                    guess[2] += deltaphi;
                    
                    if (std::abs(guess[0]) > boxsize || std::abs(guess[1]) > boxsize || std::abs(guess[2]) > boxsize)
                        return 0;
                }
                else
                    return -1;
                it++;
            }
            break;
        }
    }
    return 1;
}

int myalgorithm::getrootmultiguess(std::vector<polynomial>& poly, std::vector<double>& rhs, std::vector<double>& initialguesses, std::vector<double>& kietaphi, double boxsize, double tol, int maxit)
{
    int numguesses = initialguesses.size()/3;

    for (int i = 0; i < numguesses; i++)
    {
        kietaphi = {initialguesses[3*i+0],initialguesses[3*i+1],initialguesses[3*i+2]};
        int curstatus = getroot(poly, rhs, kietaphi, boxsize, tol, maxit);

        // Return unless Newton has not converged (if so try next initial guess):
        if (curstatus >= 0)
        {
            // if (i > 0) { std::cout << "Success at guess #" << i << std::endl; } 
            return curstatus;
        }
    }

    return -1;
}

#include "universe.h"
#include "disjointregions.h"
#include "lagrangeformfunction.h"

void myalgorithm::getreferencecoordinates(coordinategroup& coordgroup, int disjreg, std::vector<int>& elems, std::vector<double>& kietaphis)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    
    disjointregions* mydisjregs = universe::mymesh->getdisjointregions();
    elements* myelems = universe::mymesh->getelements();
    
    // Get information related to the disjoint region:
    int elemtypenum = mydisjregs->getelementtypenumber(disjreg);
    int elemdim = mydisjregs->getelementdimension(disjreg);
    int elemorder = myelems->getcurvatureorder();
    
    int rangebegin = mydisjregs->getrangebegin(disjreg), rangeend = mydisjregs->getrangeend(disjreg);
    int numelems = rangeend-rangebegin+1;

    lagrangeformfunction mylagrange(elemtypenum, elemorder, {});
    element myel(elemtypenum, elemorder);
    
    double alpha = 1.0+1.0e-8;
    if (elemorder > 1)
        alpha = 1.1;
    
    // Get the element barycenter coordinates:
    std::vector<double>* barycenters = myelems->getbarycenters(elemtypenum);
    // Get the dimensions of the box centered at the barycenter and surrounding all nodes in an element:
    std::vector<double>* boxdimensions = myelems->getboxdimensions(elemtypenum);

    // Loop on all elements in the disjoint region:
    for (int e = 0; e < numelems; e++)
    {
        double curelem = rangebegin+e;
        
        polynomials polys;
        std::vector<int> coordranking = {};
        
        double xbary = barycenters->at(3*curelem+0); double ybary = barycenters->at(3*curelem+1); double zbary = barycenters->at(3*curelem+2);
        std::vector<double> elemdist = {alpha*boxdimensions->at(3*curelem+0), alpha*boxdimensions->at(3*curelem+1), alpha*boxdimensions->at(3*curelem+2)};
        // To avoid noise related issues and to work with rotating interfaces:
        double noisedist = 0.1*(elemdist[0]+elemdist[1]+elemdist[2]);
        elemdist[0] += noisedist; elemdist[1] += noisedist; elemdist[2] += noisedist;
        
        // Loop on all candidate groups:
        coordgroup.select(xbary,ybary,zbary, *std::max_element(elemdist.begin(),elemdist.end()));
        for (int g = 0; g < coordgroup.countgroups(); g++)
        {
            std::vector<int>* curgroupindexes = coordgroup.getgroupindexes(g);
            int numcoordsingroup = curgroupindexes->size();
            if (numcoordsingroup == 0)
                continue;
            std::vector<double>* curgroupcoords = coordgroup.getgroupcoordinates(g);
            
            // Loop on all coordinates in the current group:
            for (int c = 0; c < numcoordsingroup; c++)
            {
                int curindex = curgroupindexes->at(c);
                double curx = curgroupcoords->at(3*c+0), cury = curgroupcoords->at(3*c+1), curz = curgroupcoords->at(3*c+2);

                // Only process when not yet found and when the coordinate is close enough to the element barycenter.
                if (elems[curindex] != -1 || std::abs(curx-xbary) > elemdist[0] || std::abs(cury-ybary) > elemdist[1] || std::abs(curz-zbary) > elemdist[2]) {}
                else
                {
                    // Reset initial guess:
                    std::vector<double> kietaphi = {0.0,0.0,0.0};
                    std::vector<double> rhs = {curx, cury, curz};

                    // Only create once for all coordinates the polynomials and only for the required elements:
                    if (polys.count() == 0)
                    {
                        // The coordinate polynomial used to calculate the reference coordinate must be carefully selected:
                        if (problemdimension == 3 && elemdim == 2)
                        {
                            std::vector<double> curnormal = myelems->getnormal(elemtypenum, curelem);
                            curnormal = {std::abs(curnormal[0]), std::abs(curnormal[1]), std::abs(curnormal[2])};
                            myalgorithm::stablesort(0, curnormal, coordranking);
                        }
                        else
                        {
                            myalgorithm::stablesort(0, elemdist, coordranking);
                            coordranking = {coordranking[2],coordranking[1],coordranking[0]};
                        }
                    
                        std::vector<polynomial> curpols( elemdim*(elemdim+1) );
                        for (int j = 0; j < elemdim; j++)
                        {
                            curpols[j*(elemdim+1)+0] = mylagrange.getinterpolationpolynomial(myelems->getnodecoordinates(elemtypenum, curelem, coordranking[j]));
                            for (int d = 0; d < elemdim; d++)
                                curpols[j*(elemdim+1)+1+d] = curpols[j*(elemdim+1)+0].derivative(d);
                        }
                        polys = polynomials(curpols);
                    }
                    rhs = {rhs[coordranking[0]],rhs[coordranking[1]],rhs[coordranking[2]]};
                    
                    if (myalgorithm::getroot(polys, rhs, kietaphi) == 1)
                    {
                        // Check if the (ki,eta,phi) coordinates are inside the element:
                        if (myel.isinsideelement(kietaphi[0], kietaphi[1], kietaphi[2]))
                        {
                            kietaphis[3*curindex+0] = kietaphi[0]; 
                            kietaphis[3*curindex+1] = kietaphi[1]; 
                            kietaphis[3*curindex+2] = kietaphi[2];
                            elems[curindex] = curelem;
                        }
                    }
                }
            }
        }
    }
}

std::vector<std::vector<double>> myalgorithm::splitvector(std::vector<double>& tosplit, int blocklen)
{
    int numdata = tosplit.size()/blocklen;

    std::vector<std::vector<double>> output(blocklen, std::vector<double>(numdata,0));
    
    for (int i = 0; i < numdata; i++)
    {
        for (int j = 0; j < blocklen; j++)
            output[j][i] = tosplit[i*blocklen+j];
    }
    
    return output;
}

void myalgorithm::splitvector(std::vector<int>& vec, std::vector<bool>& select, std::vector<int>& falses, std::vector<int>& trues)
{
    int numtrue = 0;
    for (int i = 0; i < select.size(); i++)
    {
        if (select[i])
            numtrue++;
    }
    
    trues = std::vector<int>(numtrue);
    falses = std::vector<int>(select.size()-numtrue);
    
    int indt = 0, indf = 0;
    for (int i = 0; i < select.size(); i++)
    {
        if (select[i])
        {
            trues[indt] = vec[i];
            indt++;
        }
        else
        {
            falses[indf] = vec[i];
            indf++;
        }
    }
}

std::vector<double> myalgorithm::normblocks(std::vector<double>& tonorm, int blocklen)
{
    int numblocks = tonorm.size()/blocklen;

    std::vector<double> output = tonorm;
    
    for (int i = 0; i < numblocks; i++)
    {
        double curnorm = 0;
        for (int j = 0; j < blocklen; j++)
            curnorm += tonorm[blocklen*i+j]*tonorm[blocklen*i+j];
        curnorm = std::sqrt(curnorm);
        
        if (curnorm > 0)
        {
            for (int j = 0; j < blocklen; j++)
                output[blocklen*i+j] = tonorm[blocklen*i+j]/curnorm;
        }
    }
    
    return output;
}

int myalgorithm::findinterval(double val, std::vector<double>& tics)
{
    int numintervals = tics.size()-1;
    // In first interval?
    if (val <= tics[1])
        return 0;
    // In last interval?
    if (val >= tics[numintervals-1])
        return numintervals-1;
    for (int i = 1; i < numintervals-1; i++)
    {
        if (val >= tics[i] && val <= tics[i+1])
            return i;
    }
}

std::string myalgorithm::getfileextension(std::string filename)
{
    int index = -1;
    for (int i = filename.length()-1; i >= 0; i--)
    {
        if (filename[i] == '.')
        {
            index = i;
            break;
        }
    }
    if (index > -1)
        return filename.substr(index, filename.length()-index);
    else
        return "";
}

std::string myalgorithm::getfilename(std::string filename)
{
    int startindex = -2, endindex = -2;
    for (int i = filename.length()-1; i >= 0; i--)
    {
        if (endindex == -2 && filename[i] == '.')
            endindex = i-1;
        if (startindex == -2 && filename[i] == '/')
            startindex = i+1;
    }
    if (startindex == -2)
        startindex = 0;
    if (endindex == -2)
        endindex = filename.length()-1;

    if (startindex > endindex)
        return "";
    else
        return filename.substr(startindex, endindex-startindex+1);
}

std::string myalgorithm::getplurals(int count)
{
    if (count > 1)
        return "s";
    else
        return "";
}

std::vector<int> myalgorithm::getequallyspaced(int start, int space, int amount)
{
    std::vector<int> output(amount);
    for (int i = 0; i < amount; i++)
        output[i] = start + space*i;

    return output;
}

void myalgorithm::osclean(std::string& line)
{
    int siz = line.size();
    if (siz > 0 && line[siz-1] == '\r')
        line.resize(siz-1);
}

std::vector<bool> myalgorithm::inttobinary(int numbits, int num)
{
    std::vector<bool> output(numbits, false);

    int index = numbits-1;
    while (num > 0)
    {
        if (num%2 == 1)
        {
            output[index] = true;
            num--;
        }
        
        num = num/2;
        index--;
    }    

    return output;
}

int myalgorithm::binarytoint(std::vector<bool> num)
{
    int output = 0;

    int pow2 = 1;
    for (int i = num.size()-1; i >= 0; i--)
    {
        if (num[i])
            output += pow2;
        pow2 = pow2*2;
    }

    return output;
}

int myalgorithm::identifyrelations(std::vector<int> numbers)
{
    if (numbers.size() <= 1)
        return 0;

    int len = numbers.size();

    int maxpos = 0;
    int fact = 1;
    for (int i = 1; i < len; i++)
    {
        if (numbers[i] > numbers[maxpos])
            maxpos = i;
        fact *= i;
    }
    
    numbers.erase(numbers.begin()+maxpos);

    return (maxpos * fact + identifyrelations(numbers));
}

int myalgorithm::factorial(int n)
{
    int out = 1;
    for (int i = 2; i <= n; i++)
        out *= i;

    return out;
}

void myalgorithm::assignedgenumbers(std::vector<std::vector<double>>& cornercoords, std::vector<int>& edgenumbers, std::vector<bool>& isbarycenteronnode, std::vector<double> noisethreshold)
{
    std::vector<int> nn(8), ne(8);
    for (int i = 0; i < 8; i++)
    {
        element myelement(i);
        nn[i] = myelement.countnodes();
        ne[i] = myelement.countedges();
    }

    // Compute the barycenter coordinates of all nodes and edges (first the edges then the nodes):
    int numnodes = 0, numedges = 0;
    for (int i = 0; i < 8; i++)
    {
        numnodes += cornercoords[i].size()/3;
        numedges += ne[i]*cornercoords[i].size()/nn[i]/3;
    }
    std::vector<double> barys(3*numedges + 3*numnodes);

    int ce = 0, cn = 0;
    for (int i = 0; i < 8; i++)
    {
        element myelement(i);
        std::vector<int> edgenodedef = myelement.getedgesdefinitionsbasedonnodes();
    
        int num = cornercoords[i].size()/nn[i]/3;
        for (int j = 0; j < num; j++)
        {
            for (int e = 0; e < ne[i]; e++)
            {
                int na = edgenodedef[2*e+0];
                int nb = edgenodedef[2*e+1];
                
                barys[3*ce+0] = 0.5*(cornercoords[i][3*nn[i]*j+3*na+0] + cornercoords[i][3*nn[i]*j+3*nb+0]);
                barys[3*ce+1] = 0.5*(cornercoords[i][3*nn[i]*j+3*na+1] + cornercoords[i][3*nn[i]*j+3*nb+1]);
                barys[3*ce+2] = 0.5*(cornercoords[i][3*nn[i]*j+3*na+2] + cornercoords[i][3*nn[i]*j+3*nb+2]);
                
                ce++;
            }
            for (int e = 0; e < nn[i]; e++)
            {
                barys[3*numedges+3*cn+0] = cornercoords[i][3*nn[i]*j+3*e+0];
                barys[3*numedges+3*cn+1] = cornercoords[i][3*nn[i]*j+3*e+1];
                barys[3*numedges+3*cn+2] = cornercoords[i][3*nn[i]*j+3*e+2];
                
                cn++;
            }
        }
    }
    
    // Sort the barycenter coordinates:
    std::vector<double> sortedbarys(barys.size());
    std::vector<int> reorderingvector;
    stablecoordinatesort(noisethreshold, barys, reorderingvector);
    for (int i = 0; i < reorderingvector.size(); i++)
    {
        sortedbarys[3*i+0] = barys[3*reorderingvector[i]+0];
        sortedbarys[3*i+1] = barys[3*reorderingvector[i]+1];
        sortedbarys[3*i+2] = barys[3*reorderingvector[i]+2];
    }
    std::vector<int> renum(reorderingvector.size());
    for (int i = 0; i < reorderingvector.size(); i++)
        renum[reorderingvector[i]] = i;
    
    // Remove duplicated barycenters:
    std::vector<int> renumberingvector;
    int numunique = removeduplicatedcoordinates(noisethreshold, sortedbarys, renumberingvector);
    
    // Assign a unique edge number for each edge:
    edgenumbers = std::vector<int>(numedges);
    for (int i = 0; i < numedges; i++)
        edgenumbers[i] = renumberingvector[renum[i]];
    
    // Calculate which edges must be split:
    std::vector<bool> isanodeatnum(numunique,false);
    for (int i = 0; i < numnodes; i++)
        isanodeatnum[renumberingvector[renum[numedges+i]]] = true;

    isbarycenteronnode = std::vector<bool>(numedges,false);
    for (int i = 0; i < numedges; i++)
    {
        if (isanodeatnum[edgenumbers[i]])
            isbarycenteronnode[i] = true;
    }
}

