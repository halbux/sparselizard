#include "gentools.h"
#include "sl.h"
#include "universe.h"
#include "disjointregions.h"
#include "lagrangeformfunction.h"
#include "slmpi.h"

#if defined(__linux__)
#include <parallel/algorithm>
#endif

void gentools::stablecoordinatesort(std::vector<double> noisethreshold, std::vector<double>& coordinates, std::vector<int>& reorderingvector)
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

void gentools::stablecoordinatesort(std::vector<double> noisethreshold, std::vector<int>& elems, std::vector<double>& coordinates, std::vector<int>& reorderingvector)
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

int gentools::removeduplicates(std::vector<double>& coordinates, std::vector<int>& renumberingvector)
{
    int numpts = coordinates.size()/3;
    renumberingvector = std::vector<int>(numpts, -1);
    
    if (numpts == 0)
        return 0;
 
    std::vector<double> noisethreshold = universe::getrawmesh()->getnodes()->getnoisethreshold();
    double ntx = noisethreshold[0], nty = noisethreshold[1], ntz = noisethreshold[2];

    coordinategroup coordgroup(coordinates);

    int numnonduplicates = 0;
    for (int i = 0; i < numpts; i++)
    {
        // If already merged there is nothing to do:
        if (renumberingvector[i] != -1)
            continue;
        renumberingvector[i] = numnonduplicates;
    
        double curx = coordinates[3*i+0], cury = coordinates[3*i+1], curz = coordinates[3*i+2];
        
        coordgroup.select(curx,cury,curz, 0);
        do
        {
            int numcoordsingroup = coordgroup.countgroupcoordinates();
            int* curgroupindexes = coordgroup.getgroupindexes();
            double* curgroupcoords = coordgroup.getgroupcoordinates();
    
            // Loop on each point in the group:
            for (int p = 0; p < numcoordsingroup; p++)
            {
                int curindex = curgroupindexes[p];
                // If close enough to be considered identical:
                if (renumberingvector[curindex] == -1 && std::abs(curx-curgroupcoords[3*p+0]) <= ntx && std::abs(cury-curgroupcoords[3*p+1]) <= nty && std::abs(curz-curgroupcoords[3*p+2]) <= ntz)
                    renumberingvector[curindex] = numnonduplicates;
            }
        }
        while (coordgroup.next());
        
        numnonduplicates++;
    }
    
    return numnonduplicates;
}

void gentools::removeduplicates(std::vector<double>& coordinates)
{
    std::vector<int> renumberingvector;

    int numuniques = removeduplicates(coordinates, renumberingvector);

    std::vector<double> uniquecoords(3*numuniques);
    for (int i = 0; i < renumberingvector.size(); i++)
    {
        int rn = renumberingvector[i];

        uniquecoords[3*rn+0] = coordinates[3*i+0];
        uniquecoords[3*rn+1] = coordinates[3*i+1];
        uniquecoords[3*rn+2] = coordinates[3*i+2];
    }

    coordinates = uniquecoords;
}

void gentools::stablesort(std::vector<int>& tosort, std::vector<int>& reorderingvector)
{
    if (reorderingvector.size() != tosort.size())
        reorderingvector.resize(tosort.size());
    
    // Set 'reorderingvector' to [0 1 2 ...]:
    std::iota(reorderingvector.begin(), reorderingvector.end(), 0);
    // Sort 'reorderingvector' according to 'tosort':
    // The < operator is overloaded by a lambda function.
    #if defined(__linux__)
    __gnu_parallel::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
    #else
    std::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
    #endif
        { 
            if (tosort[elem1] < tosort[elem2])
                return true;
            if (tosort[elem1] > tosort[elem2])
                return false;
            // For identical entries make a COHERENT decision for a stable sorting.
            return elem1 < elem2;
        });
}

void gentools::stablesort(double noisethreshold, std::vector<double>& tosort, std::vector<int>& reorderingvector)
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

void gentools::stablesort(double noisethreshold, std::vector<double>& tosort, std::vector<int>& reorderingvector, int blocklen)
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
        
    abort(); // fix return warning
}
    
void gentools::tuple3sort(std::vector<std::tuple<int,int,double>>& tosort)
{
    // Parallel sort on Linux only for now:
    #if defined(__linux__)
    __gnu_parallel::sort(tosort.begin(), tosort.end(), sortfun);
    #else
    std::sort(tosort.begin(), tosort.end(), sortfun);
    #endif
}

void gentools::slicecoordinates(std::vector<double>& toslice, double minx, double miny, double minz, double dx, double dy, double dz, int nsx, int nsy, int nsz, std::vector<int>& ga, int* pn, double* pc)
{
    int numpoints = toslice.size()/3;
    int numgroups = nsx*nsy*nsz;
    
    double invdx = 1.0/dx;
    double invdy = 1.0/dy;
    double invdz = 1.0/dz;
    
    ga = std::vector<int>(numgroups+1, 0); // one longer
    
    // Assign each coordinate to a group and count the number of coordinates in each group:
    std::vector<int> nodeingroup(numpoints);
    for (int i = 0; i < numpoints; i++)
    {
        int curxslice = std::floor((toslice[3*i+0]-minx)*invdx);
        int curyslice = std::floor((toslice[3*i+1]-miny)*invdy);
        int curzslice = std::floor((toslice[3*i+2]-minz)*invdz);
        
        // Bring in bounds:
        curxslice = std::max(curxslice,0);
        curxslice = std::min(curxslice,nsx-1);
        curyslice = std::max(curyslice,0);
        curyslice = std::min(curyslice,nsy-1);
        curzslice = std::max(curzslice,0);
        curzslice = std::min(curzslice,nsz-1);
        
        int curgroup = curxslice*nsy*nsz + curyslice*nsz + curzslice;
        
        ga[curgroup+1]++; // +1
        nodeingroup[i] = curgroup;
    }

    // Assign the group ranges (ga[0] is zero):
    for (int g = 1; g < numgroups+1; g++)
        ga[g] = ga[g-1] + ga[g];
    
    // Populate the remaining outputs:
    std::vector<int> indexingroup(numgroups,0);
    for (int i = 0; i < numpoints; i++)
    {
        int curgroup = nodeingroup[i];
        int curpos = ga[curgroup] + indexingroup[curgroup];
        
        pn[curpos] = i;
        pc[3*curpos+0] = toslice[3*i+0];
        pc[3*curpos+1] = toslice[3*i+1];
        pc[3*curpos+2] = toslice[3*i+2];
        
        indexingroup[curgroup]++;
    }
}

std::vector<double> gentools::getcoordbounds(std::vector<double>& coords)
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

std::vector<int> gentools::unique(std::vector<int> a)
{
    std::sort(a.begin(), a.end());
    std::vector<int>::iterator it;
    it = std::unique(a.begin(), a.end());
    a.resize(std::distance(a.begin(), it));
    
    return a;
}

std::vector<int> gentools::intersect(std::vector<int> a, std::vector<int> b)
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

std::vector<int> gentools::exclude(std::vector<int> a, std::vector<int> b)
{
    // Sort both vectors:
    std::sort(a.begin(), a.end()); 
    std::sort(b.begin(), b.end());
    
    std::vector<int> output(a.size()+b.size());
    std::vector<int>::iterator it;
    
    // Exclude:
    it = std::set_difference(a.begin(), a.end(), b.begin(), b.end(), output.begin());
    output.resize(it-output.begin()); 
    
    return output;
}

void gentools::csrtoijk(int numberofrows, int* csrrows, int* ijkrows)
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

int gentools::getroot(polynomials& polys, std::vector<double>& rhs, std::vector<double>& guess, double boxsize, double tol, int maxit)
{
    int it = 0;
    double deltaki = 1.0, deltaeta = 1.0, deltaphi = 1.0;
    
    // Limit the jumps to a fraction of the box size:
    double maxjump = 0.45;

    switch (polys.count())
    {
        case 1:
        {
            while (std::abs(deltaki) > tol)
            {
                if (it < maxit)
                {
                    std::vector<double> evaled;
                    polys.evalatsingle(guess, 1, evaled);
                
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
        case 2:
        {
            while (std::abs(deltaki) > tol || std::abs(deltaeta) > tol)
            {
                if (it < maxit)
                {
                    std::vector<double> evaled;
                    polys.evalatsingle(guess, 2, evaled);
                    
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
        case 3:
        {
            while (std::abs(deltaki) > tol || std::abs(deltaeta) > tol || std::abs(deltaphi) > tol)
            {
                if (it < maxit)
                {
                    std::vector<double> evaled;
                    polys.evalatsingle(guess, 3, evaled);
                    
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

void gentools::getreferencecoordinates(coordinategroup& coordgroup, int disjreg, std::vector<int>& elems, std::vector<double>& kietaphis)
{
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    
    disjointregions* mydisjregs = universe::getrawmesh()->getdisjointregions();
    elements* myelems = universe::getrawmesh()->getelements();
    
    // Get information related to the disjoint region:
    int elemtypenum = mydisjregs->getelementtypenumber(disjreg);
    int elemdim = mydisjregs->getelementdimension(disjreg);
    int elemorder = myelems->getcurvatureorder();
    
    int rangebegin = mydisjregs->getrangebegin(disjreg), rangeend = mydisjregs->getrangeend(disjreg);
    int numelems = rangeend-rangebegin+1;

    polynomials polys(lagrangeformfunction(elemtypenum,elemorder,{}).getformfunctionpolynomials());
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
        
        polynomials syspolys;
        std::vector<int> coordranking = {};
        
        double xbary = barycenters->at(3*curelem+0); double ybary = barycenters->at(3*curelem+1); double zbary = barycenters->at(3*curelem+2);
        std::vector<double> elemdist = {alpha*boxdimensions->at(3*curelem+0), alpha*boxdimensions->at(3*curelem+1), alpha*boxdimensions->at(3*curelem+2)};
        // To avoid noise related issues and to work with rotating interfaces:
        double noisedist = 0.1*(elemdist[0]+elemdist[1]+elemdist[2]);
        elemdist[0] += noisedist; elemdist[1] += noisedist; elemdist[2] += noisedist;
        
        // Loop on all candidate groups:
        coordgroup.select(xbary,ybary,zbary, *std::max_element(elemdist.begin(),elemdist.end()));
        do
        {
            int numcoordsingroup = coordgroup.countgroupcoordinates();
            int* curgroupindexes = coordgroup.getgroupindexes();
            double* curgroupcoords = coordgroup.getgroupcoordinates();
            
            // Loop on all coordinates in the current group:
            for (int c = 0; c < numcoordsingroup; c++)
            {
                int curindex = curgroupindexes[c];
                double curx = curgroupcoords[3*c+0], cury = curgroupcoords[3*c+1], curz = curgroupcoords[3*c+2];

                // Only process when not yet found and when the coordinate is close enough to the element barycenter.
                if (elems[curindex] != -1 || std::abs(curx-xbary) > elemdist[0] || std::abs(cury-ybary) > elemdist[1] || std::abs(curz-zbary) > elemdist[2]) {}
                else
                {
                    // Reset initial guess:
                    std::vector<double> kietaphi = {0.0,0.0,0.0};
                    std::vector<double> rhs = {curx, cury, curz};

                    // Only create once for all coordinates the polynomials and only for the required elements:
                    if (syspolys.count() == 0)
                    {
                        // The coordinate polynomial used to calculate the reference coordinate must be carefully selected:
                        if (problemdimension == 3 && elemdim == 2)
                        {
                            std::vector<double> curnormal = myelems->getnormal(elemtypenum, curelem);
                            curnormal = {std::abs(curnormal[0]), std::abs(curnormal[1]), std::abs(curnormal[2])};
                            stablesort(0, curnormal, coordranking);
                        }
                        else
                        {
                            stablesort(0, elemdist, coordranking);
                            coordranking = {coordranking[2],coordranking[1],coordranking[0]};
                        }
                        
                        std::vector<int> trimmedcr = coordranking;
                        trimmedcr.resize(elemdim);
                    
                        std::vector<double> curcoords = myelems->getnodecoordinates(elemtypenum, curelem);
                        
                        std::vector<double> xyz = gentools::separate(curcoords, 3, trimmedcr);
                        syspolys = polys.sum(xyz);
                    }
                    rhs = {rhs[coordranking[0]],rhs[coordranking[1]],rhs[coordranking[2]]};
                    
                    if (getroot(syspolys, rhs, kietaphi) == 1)
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
        while (coordgroup.next());
    }
}

std::vector<std::vector<double>> gentools::splitvector(std::vector<double>& tosplit, int blocklen)
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

void gentools::splitvector(std::vector<int>& vec, std::vector<bool>& select, std::vector<int>& falses, std::vector<int>& trues)
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

void gentools::select(std::vector<int>& vals, std::vector<int>& selectedindexes, std::vector<int>& selected)
{
    int numselected = selectedindexes.size();
    selected.resize(numselected);

    for (int i = 0; i < numselected; i++)
        selected[i] = vals[selectedindexes[i]];
}

void gentools::select(std::vector<bool>& vals, indexmat selectedindexes, std::vector<bool>& selected)
{
    int numselected = selectedindexes.count();
    selected.resize(numselected);

    int* selptr = selectedindexes.getvalues();

    for (int i = 0; i < numselected; i++)
        selected[i] = vals[selptr[i]];
}

bool gentools::isflipped(std::vector<int>& a, std::vector<int>& b)
{
    int num = a.size();
    if (num == 1)
        return false;
    if (num == 2)
        return (a[0] != b[0]);

    // Find the corresponding node in b of a[0]:
    int ba0 = -1;
    for (int i = 0; i < num; i++)
    {
        if (b[i] == a[0])
        {
            ba0 = i;
            break;
        }
    }
    
    // Decide with position of neighbour:
    if (b[(ba0+1)%num] == a[1])
        return false;
    else
        return true;
}

std::vector<double> gentools::normblocks(std::vector<double>& tonorm, int blocklen)
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

int gentools::findinterval(double val, std::vector<double>& tics)
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
    
    abort(); // fix return warning
}

std::vector<double> gentools::getintervaltics(double minval, double maxval, int numintervals)
{
    double range = maxval - minval;
    double step = range/((double)numintervals);

    std::vector<double> output(numintervals+1, minval);
    for (int i = 1; i < numintervals+1; i++)
        output[i] = output[i-1] + step;
    
    return output;
}

std::string gentools::getfileextension(std::string filename)
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

std::string gentools::getfilename(std::string filename)
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

std::string gentools::getplurals(int count)
{
    if (count > 1)
        return "s";
    else
        return "";
}

std::vector<int> gentools::getequallyspaced(int start, int space, int amount)
{
    std::vector<int> output(amount);
    for (int i = 0; i < amount; i++)
        output[i] = start + space*i;

    return output;
}

std::vector<double> gentools::duplicate(std::vector<double> invec, int n)
{
    int len = invec.size();
    std::vector<double> out(len*n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < len; j++)
            out[len*i+j] = invec[j];
    }
    return out;
}

void gentools::osclean(std::string& line)
{
    int siz = line.size();
    if (siz > 0 && line[siz-1] == '\r')
        line.resize(siz-1);
}

std::vector<bool> gentools::inttobinary(int numbits, int num)
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

int gentools::binarytoint(std::vector<bool> num)
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

double gentools::exactinttodouble(long long int num)
{
    // All integers of absolute value < 2^53 can be exactly represented in IEEE double format:
    if (std::abs(num) < 9007199254740992)
        return (double)num;
    else
    {
        std::cout << "Error in 'gentools' namespace: integer " << num << " is too large to be exactly represented in double format" << std::endl;
        abort();
    }
}

int gentools::identifyrelations(std::vector<int> numbers)
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

int gentools::factorial(int n)
{
    int out = 1;
    for (int i = 2; i <= n; i++)
        out *= i;

    return out;
}

void gentools::assignedgenumbers(std::vector<bool>& isownelem, std::vector<std::vector<double>>& cornercoords, std::vector<int>& edgenumbers, std::vector<bool>& isbarycenteronnode)
{
    int numranks = slmpi::count();
    int rank = slmpi::getrank();
    
    std::shared_ptr<dtracker> dt = universe::getrawmesh()->getdtracker();
    int numneighbours = dt->countneighbours();
    std::vector<int> neighbours = dt->getneighbours();
    
    std::vector<int> nn(8), ne(8);
    for (int i = 0; i < 8; i++)
    {
        element myelement(i);
        nn[i] = myelement.countnodes();
        ne[i] = myelement.countedges();
    }

    // Compute the barycenter coordinates of all edges:
    int numedges = 0;
    for (int i = 0; i < 8; i++)
        numedges += ne[i]*cornercoords[i].size()/nn[i]/3;
        
    int numownedges = 0;
    std::vector<bool> isownedge(numedges, false);
    std::vector<double> barys(3*numedges);

    int ce = 0, cc = 0;
    for (int i = 0; i < 8; i++)
    {
        element myelement(i);
        std::vector<int> edgenodedef = myelement.getedgesdefinitionsbasedonnodes();
    
        int num = cornercoords[i].size()/nn[i]/3;
        for (int j = 0; j < num; j++)
        {
            for (int e = 0; e < ne[i]; e++)
            {
                // FOR PYRAMIDS THIS IS INCORRECT SINCE FULLSPLIT PYRAMIDS CREATE TWO ELEMENT TYPES.
                // THUS THE [i][j] ORDERING DOES NOT FOLLOW LEAVES FOR INCREASING J.
                if (isownelem[cc])
                {
                    numownedges++;
                    isownedge[ce] = true;
                }
                
                int na = edgenodedef[2*e+0];
                int nb = edgenodedef[2*e+1];
                
                barys[3*ce+0] = 0.5*(cornercoords[i][3*nn[i]*j+3*na+0] + cornercoords[i][3*nn[i]*j+3*nb+0]);
                barys[3*ce+1] = 0.5*(cornercoords[i][3*nn[i]*j+3*na+1] + cornercoords[i][3*nn[i]*j+3*nb+1]);
                barys[3*ce+2] = 0.5*(cornercoords[i][3*nn[i]*j+3*na+2] + cornercoords[i][3*nn[i]*j+3*nb+2]);
                
                ce++;
            }
            
            cc++;
        }
    }
    
    // Append the corner nodes to the barycenters:
    std::vector<double> catcornercoords;
    concatenate(cornercoords, catcornercoords);
    removeduplicates(catcornercoords);
    std::vector<int> neighboursnumownedges = appendneighbourvalues(barys, catcornercoords, numownedges);
    
    // Find duplicates:
    std::vector<int> renumberingvector;
    int numunique = removeduplicates(barys, renumberingvector);
    
    // Calculate which edges must be split:
    std::vector<bool> isanodeatnum(numunique, false);
    for (int i = 0; i < barys.size()/3-numedges; i++)
        isanodeatnum[renumberingvector[numedges+i]] = true;

    isbarycenteronnode = std::vector<bool>(numedges, false);
    for (int i = 0; i < numedges; i++)
    {
        if (isanodeatnum[renumberingvector[i]])
            isbarycenteronnode[i] = true;
    }
    
    // Assign a unique edge number for each edge:
    edgenumbers = std::vector<int>(numedges);
    for (int i = 0; i < numedges; i++)
        edgenumbers[i] = renumberingvector[i];

    // Make the edges numbers continuous:
    std::vector<int> edgerenum;
    int numuniqueedges = squeeze(edgenumbers, numunique, edgerenum);
    for (int i = 0; i < numedges; i++)
        edgenumbers[i] = edgerenum[edgenumbers[i]];
    
    // Harmonize the edges numbers across neighbour ranks:
    if (dt->isdefined() == false || numranks == 1)
        return;
        
    std::vector<int> ownedgesindexes;
    find(isownedge, numownedges, ownedgesindexes);
        
    barys.resize(3*numedges);
    
    std::vector<double> ownbarys(3*numownedges);
    selectcoordinates(isownedge, barys, ownbarys.data());
    
    // Get the edges count on all ranks to create unique global edge numbers:
    std::vector<int> fragment = {numuniqueedges};
    std::vector<int> allnumuniqueedges;
    slmpi::allgather(fragment, allnumuniqueedges);
    
    // USE LONG LONG INT TO ALLOW MORE EDGES IN THE MESH.
    std::vector<int> edgenumshift(numranks, 0);
    for (int i = 1; i < numranks; i++)
        edgenumshift[i] = edgenumshift[i-1] + allnumuniqueedges[i-1];
    
    for (int i = 0; i < numedges; i++)
        edgenumbers[i] += edgenumshift[rank];
    
    // Extract the own edges numbers:
    std::vector<int> ownedgenumbers;
    select(edgenumbers, ownedgesindexes, ownedgenumbers);
    
    std::vector<std::vector<double>> ownbaryssends(numneighbours, ownbarys);
    std::vector<std::vector<double>> ownbarysrecvs(numneighbours);
    for (int n = 0; n < numneighbours; n++)
        ownbarysrecvs[n].resize(3*neighboursnumownedges[n]);
        
    slmpi::exchange(neighbours, ownbaryssends, ownbarysrecvs);
    
    // Harmonize the no-overlap interface:
    std::vector<std::vector<int>> ownedgenumberssends(numneighbours, ownedgenumbers);
    std::vector<std::vector<int>> ownedgenumbersrecvs(numneighbours);
    for (int n = 0; n < numneighbours; n++)
        ownedgenumbersrecvs[n].resize(neighboursnumownedges[n]);
        
    slmpi::exchange(neighbours, ownedgenumberssends, ownedgenumbersrecvs);

    // Reused below:
    std::vector<std::vector<int>> posfound(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        gentools::findcoordinates(ownbarysrecvs[n], barys, posfound[n]);

        for (int i = 0; i < posfound[n].size(); i++)
        {
            if (posfound[n][i] != -1 && isownedge[i])
                edgenumbers[i] = std::min(edgenumbers[i], ownedgenumbersrecvs[n][posfound[n][i]]); // min rank decides
        }
    }

    // Harmonize the outer overlap interface:
    if (dt->isoverlap() == false)
        return;

    // Update the own edges numbers:
    select(edgenumbers, ownedgesindexes, ownedgenumbers);

    // Send the no-overlap harmonized edge numbers:
    ownedgenumberssends = std::vector<std::vector<int>>(numneighbours, ownedgenumbers);
    slmpi::exchange(neighbours, ownedgenumberssends, ownedgenumbersrecvs);

    for (int n = 0; n < numneighbours; n++)
    {
        for (int i = 0; i < posfound[n].size(); i++)
        {
            if (posfound[n][i] != -1)
                edgenumbers[i] = ownedgenumbersrecvs[n][posfound[n][i]];
        }
    }
}

std::vector<double> gentools::separate(std::vector<double>& v, int blocklen, std::vector<int> sel)
{
    int numblocks = v.size()/blocklen;
    
    std::vector<double> output(sel.size() * numblocks);
    
    for (int s = 0; s < sel.size(); s++)
    {
        for (int b = 0; b < numblocks; b++)
            output[s*numblocks+b] = v[b*blocklen + sel[s]];
    }
    
    return output;
}

std::vector<int> gentools::chainrenumbering(std::vector<int>& originalrenum, std::vector<int>& newrenum)
{
    std::vector<int> output(originalrenum.size());

    for (int i = 0; i < originalrenum.size(); i++)
        output[i] = newrenum[originalrenum[i]];

    return output;
}

std::vector<int> gentools::invertrenumbering(std::vector<int>& renum)
{
    std::vector<int> output(renum.size());

    for (int i = 0; i < renum.size(); i++)
        output[renum[i]] = i;

    return output;
}

std::vector<int> gentools::getreordering(std::vector<int>& renum)
{
    std::vector<int> out(renum.size());
    
    for (int i = 0; i < renum.size(); i++)
        out[renum[i]] = i;
    
    return out;
}

void gentools::reorder(std::vector<int>& inad, std::vector<double>& indat, std::vector<int>& renumbering, std::vector<int>& outad, std::vector<double>& outdat)
{
    if (indat.size() == 0)
    {
        outad = {};
        outdat = {};
        return;
    }

    int num = inad.size()-1;
    
    std::vector<int> outcnt(num,0);
    for (int i = 0; i < num; i++)
        outcnt[renumbering[i]] = inad[i+1]-inad[i];
        
    outad = std::vector<int>(num+1,0);
    for (int i = 1; i < num+1; i++)
        outad[i] = outad[i-1] + outcnt[i-1];  

    outdat = std::vector<double>(indat.size());
    for (int i = 0; i < num; i++)
    {
        int inpos = inad[i];
        int outpos = outad[renumbering[i]];
        for (int j = 0; j < inad[i+1]-inad[i]; j++)
            outdat[outpos+j] = indat[inpos+j];
    }
}

void gentools::toaddressdata(std::vector<int>& elems, std::vector<double>& refcoords, std::vector<int> totalnumelems, std::vector<std::vector<int>>& ads, std::vector<std::vector<double>>& rcs, std::vector<int>& indexinrcsoforigin)
{
    int ne = elems.size()/2;

    std::vector<int> countintype(8,0);
    std::vector<std::vector<int>> cnt(8);
    for (int i = 0; i < 8; i++)
        cnt[i] = std::vector<int>(totalnumelems[i],0);
        
    for (int i = 0; i < ne; i++)
    {
        int typ = elems[2*i+0];
        int num = elems[2*i+1];
        
        countintype[typ]++;
        cnt[typ][num]++;
    }
    
    // Remove empty types:
    std::vector<std::vector<int>> index(8);
    for (int i = 0; i < 8; i++)
    {
        if (countintype[i] == 0)
            cnt[i] = {};
        index[i] = std::vector<int>(cnt[i].size(),0);
    }
    
    ads = std::vector<std::vector<int>>(8, std::vector<int>(0));
    rcs = std::vector<std::vector<double>>(8, std::vector<double>(0));
    
    for (int i = 0; i < 8; i++)
    {
        if (cnt[i].size() == 0)
            continue;
    
        ads[i] = std::vector<int>(cnt[i].size()+1,0);
        for (int j = 1; j < ads[i].size(); j++)
            ads[i][j] = ads[i][j-1] + 3*cnt[i][j-1];
            
        rcs[i] = std::vector<double>(3*countintype[i]);
    }
    
    indexinrcsoforigin = std::vector<int>(ne);
    for (int i = 0; i < ne; i++)
    {
        int typ = elems[2*i+0];
        int num = elems[2*i+1];
        int pos = ads[typ][num] + index[typ][num];
        
        rcs[typ][pos+0] = refcoords[3*i+0];
        rcs[typ][pos+1] = refcoords[3*i+1];
        rcs[typ][pos+2] = refcoords[3*i+2];
        
        indexinrcsoforigin[i] = pos/3;
        
        index[typ][num] += 3;
    }
}

std::vector<int> gentools::concatenate(std::vector<std::vector<int>> tocat)
{
    // Get the total size:
    int len = 0;
    for (int i = 0; i < tocat.size(); i++)
        len += tocat[i].size();
        
    std::vector<int> output(len);
    int index = 0;
    for (int i = 0; i < tocat.size(); i++)
    {
        for (int j = 0; j < tocat[i].size(); j++)
        {
            output[index] = tocat[i][j];
            index++;
        }
    }
    
    return output;
}

void gentools::concatenate(std::vector<std::vector<double>>& tocat, std::vector<double>& cated)
{
    // Get the total size:
    int len = 0;
    for (int i = 0; i < tocat.size(); i++)
        len += tocat[i].size();

    cated.resize(len);

    int index = 0;
    for (int i = 0; i < tocat.size(); i++)
    {
        for (int j = 0; j < tocat[i].size(); j++)
            cated[index+j] = tocat[i][j];

        index += tocat[i].size();
    }
}

int gentools::inequalitytoint(int a, int b)
{
    if (a < b)
        return -1;
    else
    {
        if (a == b)
            return 0;
        else
            return 1;
    }
}

void gentools::normvector(std::vector<double>& tonorm)
{
    // Compute the norm:
    double nrm = 0.0;
    for (int i = 0; i < tonorm.size(); i++)
        nrm += tonorm[i]*tonorm[i];
    nrm = std::sqrt(nrm);
    double invnrm = 1.0/nrm;

    for (int i = 0; i < tonorm.size(); i++)
        tonorm[i] = invnrm * tonorm[i];
}

void gentools::solveuppertriangular(int len, double* U, double* b, double* x)
{
    // Init x to all zero:
    for (int i = 0; i < len; i++)
        x[i] = 0;
 
    int index = ((1+len)*len)/2 - 1;
    
    for (int c = len-1; c >= 0; c--)
    {   
        x[c] = (b[c]-x[c]) / U[index];
        
        index--;
        
        for (int r = c-1; r >= 0; r--)
        {
            x[r] += U[index] * x[c];
            
            index--;
        }
    }
}

void gentools::givensrotation(double a, double b, double& c, double& s, double& r)
{
    // As in www.netlib.org/eispack/comqr.f
    r = std::sqrt(a*a + b*b);
    c = a/r;
    s = b/r;
}

void gentools::applygivensrotation(double* h, std::vector<double>& cs, std::vector<double>& sn, int k)
{
    for (int i = 0; i < k; i++)
    {
        double temp = cs[i] * h[i] + sn[i] * h[i+1];
        h[i+1] = -sn[i] * h[i] + cs[i] * h[i+1];
        h[i] = temp;
    }

    // Update the next sin and cos values for the rotation:
    double r;
    givensrotation(h[k], h[k+1], cs[k], sn[k], r);

    // Eliminate h[k+1]:
    h[k] = r;
    h[k+1] = 0.0;
}

std::vector<double> gentools::arnoldi(densemat (*mymatmult)(densemat), densemat Q, int k)
{   
    // One Krylov vector on each row:
    int n = Q.countcolumns();
    double* Qptr = Q.getvalues();
    
    // Krylov vector fragment on each rank:
    densemat q = mymatmult(Q.extractrows(k,k).getresized(n,1));
    double* qptr = q.getvalues();
    
    if (q.countrows() != n || q.countcolumns() != 1)
    {
        std::cout << "Error in 'gentools' namespace: in function arnoldi the matrix product function call returned a densemat of wrong size on rank " << slmpi::getrank() << std::endl;
        abort();
    }
    
    // Standard Gramm-Schmidt orthogonalization:
    densemat h = Q.getresized(k+1,n).multiply(q);
    h = h.getresized(1,k+1);
    
    slmpi::sum(k+1, h.getvalues()); // reduce on all ranks
    
    densemat Qh = h.multiply(Q.getresized(k+1,n));
    double* Qhptr = Qh.getvalues();
    
    // Subtract Qh and compute the norm of q:
    double normqfrag = 0.0;
    for (int i = 0; i < n; i++)
    {
        qptr[i] -= Qhptr[i];
        normqfrag += qptr[i]*qptr[i];
    }
        
    slmpi::sum(1, &normqfrag); // reduce on all ranks
    
    double normq = std::sqrt(normqfrag);

    // Norm the Krylov vector and place it in Q:
    double invnormq = 1.0/normq;
    for (int i = 0; i < n; i++)
        Qptr[(k+1)*n+i] = invnormq * qptr[i];
    
    std::vector<double> hvec;
    h.getvalues(hvec);
    hvec.push_back(normq);
    
    return hvec;
}

int gentools::squeeze(std::vector<int>& nums, int maxval, std::vector<int>& renumbering)
{
    std::vector<bool> isactive(maxval+1, false);
    for (int i = 0; i < nums.size(); i++)
        isactive[nums[i]] = true;

    int index = 0;
    renumbering = std::vector<int>(maxval+1, -1);
    for (int i = 0; i < maxval+1; i++)
    {
        if (isactive[i])
        {
            renumbering[i] = index;
            index++;
        }
    }

    return index;
}

void gentools::find(std::vector<bool>& invec, int numtrue, std::vector<int>& trueindexes)
{
    if (numtrue == -1)
    {
        numtrue = 0;
        for (int i = 0; i < invec.size(); i++)
        {
            if (invec[i])
                numtrue++;
        }
    }
    
    trueindexes = std::vector<int>(numtrue);
    
    int ind = 0;
    for (int i = 0; i < invec.size(); i++)
    {
        if (invec[i])
        {
            trueindexes[ind] = i;
            ind++;
        }
    }
}

int gentools::findcoordinates(std::vector<double>& targetcoords, std::vector<double>& tofindintarget, std::vector<int>& posfound)
{
    int nct = targetcoords.size()/3;
    int ncf = tofindintarget.size()/3;

    // Place all in one vector:
    std::vector<double> allcoords(targetcoords.size()+tofindintarget.size());
    for (int i = 0; i < targetcoords.size(); i++)
        allcoords[i] = targetcoords[i];
    for (int i = 0; i < tofindintarget.size(); i++)
        allcoords[targetcoords.size()+i] = tofindintarget[i];

    std::vector<int> renumberingvector;
    gentools::removeduplicates(allcoords, renumberingvector);

    // This is needed for general renumbering vectors (also allows duplicates in the target):
    std::vector<int> targetpositions(nct+ncf, -1);
    for (int i = 0; i < nct; i++)
        targetpositions[renumberingvector[i]] = i;

    int numfound = 0;
    posfound = std::vector<int>(ncf, -1);
    for (int i = 0; i < ncf; i++)
    {
        int renumed = renumberingvector[nct+i];
        if (targetpositions[renumed] != -1)
        {
            posfound[i] = targetpositions[renumed];
            numfound++;
        }
    }

    return numfound;
}

void gentools::selectcoordinates(std::vector<bool>& selection, std::vector<double>& coords, double* selectedcoords)
{
    int index = 0;
    for (int i = 0; i < selection.size(); i++)
    {
        if (selection[i])
        {
            selectedcoords[3*index+0] = coords[3*i+0];
            selectedcoords[3*index+1] = coords[3*i+1];
            selectedcoords[3*index+2] = coords[3*i+2];
            
            index++;
        }
    }
}

void gentools::pickcandidates(int numbertopick, std::vector<double>& candidatecoordinates, std::vector<double>& picked)
{
    picked = {};
    
    int numcandidates = candidatecoordinates.size()/3;
    
    if (numbertopick <= 0 || numcandidates == 0)
        return;

    picked = std::vector<double>(3*numbertopick);

    // Pick should include all candidates if the number to pick is larger than the number of candidates: 
    int spacing = ceildiv(numcandidates, numbertopick);

    for (int i = 0; i < numbertopick; i++)
    {
        int p = i*spacing;
        if (p >= numcandidates)
            p = numcandidates-1;
    
        picked[3*i+0] = candidatecoordinates[3*p+0];
        picked[3*i+1] = candidatecoordinates[3*p+1];
        picked[3*i+2] = candidatecoordinates[3*p+2];
    }
}

void gentools::split(std::vector< std::vector<std::vector<double>> >& data, std::vector<double>& dataa, std::vector<double>& datab, std::vector<std::vector<int>>& sizesa, std::vector<std::vector<int>>& sizesb)
{
    int sizea = 0;
    int sizeb = 0;

    sizesa.resize(data.size());
    sizesb.resize(data.size());

    for (int i = 0; i < data.size(); i++)
    {
        sizesa[i] = std::vector<int>(data[i].size(), 0);
        sizesb[i] = std::vector<int>(data[i].size(), 0);

        for (int j = 0; j < data[i].size(); j++)
        {
            if (data[i][j].size() >= 2)
            {
                int csa = (int)data[i][j][0];
                int csb = (int)data[i][j][1];

                sizesa[i][j] = csa;
                sizesb[i][j] = csb;

                sizea += csa;
                sizeb += csb;
            }
        }
    }

    dataa = std::vector<double>(sizea);
    datab = std::vector<double>(sizeb);

    int indexa = 0;
    int indexb = 0;

    for (int i = 0; i < data.size(); i++)
    {
        for (int j = 0; j < data[i].size(); j++)
        {
            if (data[i][j].size() >= 2)
            {
                int csa = (int)data[i][j][0];
                int csb = (int)data[i][j][1];

                for (int k = 0; k < csa; k++)
                    dataa[indexa+k] = data[i][j][2+k];

                for (int k = 0; k < csb; k++)
                    datab[indexb+k] = data[i][j][2+csa + k];

                indexa += csa;
                indexb += csb;
            }
        }
    }
}

void gentools::pack(std::vector<int> tags, std::vector<std::vector<double>>& topack, std::vector<std::vector<double>>& packed)
{
    int numtags = tags.size();

    // Get the packed sizes:
    std::vector<int> packedsizes(numtags, 0);
    for (int i = 0; i < numtags; i++)
    {
        for (int j = i+1; j < numtags; j++)
        {
            int ci = i*numtags+j;
            int cs = topack[ci].size();
         
            if (cs == 0)
                continue;
                
            packedsizes[i] += 2 + cs;
            packedsizes[j] += 2 + cs;
        }
    }

    // Preallocate:
    packed = std::vector<std::vector<double>>(numtags, std::vector<double>(0));
    for (int i = 0; i < numtags; i++)
        packed[i] = std::vector<double>(packedsizes[i]);

    // Populate:
    std::vector<int> indexes(numtags, 0);
  
    for (int i = 0; i < numtags; i++)
    {
        for (int j = i+1; j < numtags; j++)
        {
            int ci = i*numtags+j;
            int cs = topack[ci].size();
            
            if (cs == 0)
                continue;

            packed[i][indexes[i]+0] = exactinttodouble(tags[j]);
            packed[i][indexes[i]+1] = exactinttodouble(cs);

            packed[j][indexes[j]+0] = exactinttodouble(tags[i]);
            packed[j][indexes[j]+1] = exactinttodouble(cs);

            for (int k = 0; k < cs; k++)
            {
                double val = topack[ci][k];
                packed[i][indexes[i]+2+k] = val;
                packed[j][indexes[j]+2+k] = val;
            }

            indexes[i] += 2+cs;
            indexes[j] += 2+cs;
        }
    }
}

std::vector<int> gentools::unpack(std::vector<double>& packed, std::vector<std::vector<double>>& unpacked)
{
    // Get the number of datasets:
    int numdatasets = 0;
    
    int index = 0;
    while (index < packed.size())
    {
        index += 2 + (int)packed[index+1];
        numdatasets++;
    }

    std::vector<int> output(numdatasets);
    unpacked = std::vector<std::vector<double>>(numdatasets, std::vector<double>(0));

    index = 0;
    int cds = 0;
    while (index < packed.size())
    {
        output[cds] = (int)packed[index];
        int len = (int)packed[index+1];
        
        index += 2;
        
        unpacked[cds] = std::vector<double>(len);
        double* data = unpacked[cds].data();
        
        for (int i = 0; i < len; i++)
            data[i] = packed[index+i];

        index += len;
        cds++;
    }

    return output;
}

std::vector<int> gentools::extract(std::vector<int>& data, int period, int shift)
{
    int numperiods = data.size()/period;
    
    std::vector<int> output(numperiods);
    
    int index = 0;
    for (int i = 0; i < numperiods; i++)
    {
        for (int j = 0; j < period; j++)
        {
            if (j == shift)
                output[i] = data[i*period+j];
            else
            {
                data[index] = data[i*period+j];
                index++;
            }
        }
    }
    data.resize(index);

    return output;
}

std::vector<double> gentools::extract(std::vector<double>& data, int period, int shift)
{
    int numperiods = data.size()/period;
    
    std::vector<double> output(numperiods);
    
    int index = 0;
    for (int i = 0; i < numperiods; i++)
    {
        for (int j = 0; j < period; j++)
        {
            if (j == shift)
                output[i] = data[i*period+j];
            else
            {
                data[index] = data[i*period+j];
                index++;
            }
        }
    }
    data.resize(index);

    return output;
}

int gentools::ceildiv(int a, int b)
{
    int r = a%b;
    if (r == 0)
        return a/b;
    else
        return (a-r)/b+1;
}

int gentools::getpackedsize(int numbits)
{
    int numbitsinint = 8 * sizeof(int); // to be os independent
    numbitsinint--; // sign bit is not used to avoid issues on 0 int value
    
    return ceildiv(numbits, numbitsinint);
}

void gentools::pack(std::vector<bool>& topack, std::vector<int>& packed)
{   
    int numbitsinint = 8 * sizeof(int); // to be os independent
    numbitsinint--; // sign bit is not used to avoid issues on 0 int value
 
    int numbits = topack.size();
    int numints = ceildiv(numbits, numbitsinint);
    
    packed = std::vector<int>(numints);
    
    std::vector<int> powersof2(numbitsinint, 1);
    for (int i = 1; i < numbitsinint; i++)
        powersof2[i] = 2*powersof2[i-1];
    
    for (int i = 0; i < numints; i++)
    {
        int val = 0;
        for (int b = 0; b < std::min(numbitsinint, numbits-i*numbitsinint); b++)
        {
            if (topack[i*numbitsinint+b])
                val += powersof2[b];
        }
        packed[i] = val;
    }
}

void gentools::unpack(int orignumbools, std::vector<int>& packed, std::vector<bool>& unpacked)
{
    unpacked = std::vector<bool>(orignumbools, false);

    int numbitsinint = 8 * sizeof(int); // to be os independent
    numbitsinint--; // sign bit is not used to avoid issues on 0 int value
 
    int numints = packed.size();
 
    for (int i = 0; i < numints; i++)
    {
        int val = packed[i];
        for (int b = 0; b < numbitsinint; b++)
        {
            if (val > 0)
            {
                if (val%2 == 1)
                {
                    unpacked[i*numbitsinint+b] = true;
                    val--;
                }
                val = val/2;
            }
        }
    }
}

void gentools::split(std::vector<int>& vals, int val, std::vector<int>& lowervals, std::vector<int>& highervals)
{
    int nl = 0, nh = 0;
    for (int i = 0; i < vals.size(); i++)
    {
        if (vals[i] < val)
            nl++;
        else
            nh++;
    }
    
    lowervals = std::vector<int>(nl);
    highervals = std::vector<int>(nh);
    
    int il = 0, ih = 0;
    for (int i = 0; i < vals.size(); i++)
    {
        if (vals[i] < val)
        {
            lowervals[il] = vals[i];
            il++;
        }
        else
        {
            highervals[ih] = vals[i];
            ih++;
        }
    }
}

int gentools::counttrue(std::vector<bool>& tocount)
{
    int cnt = 0;
    for (int i = 0; i < tocount.size(); i++)
    {
        if (tocount[i])
            cnt++;
    }
    return cnt;
}

void gentools::compresszeros(std::vector<int>& tocompress)
{
    int len = tocompress.size();
    
    int i = 0, ci = 0;
    while (i < len)
    {
        if (tocompress[i] > 0)
        {
            tocompress[ci] = tocompress[i];
            i++; ci++;
        }
        else
        {
            int cntconsec = 0;
            while (i < len && tocompress[i] == 0)
            {
                cntconsec++;
                i++;
            }
            tocompress[ci] = -cntconsec;
            ci++;
        }
    }
    
    tocompress.resize(ci+1);
    tocompress[ci] = len;
}

void gentools::decompresszeros(std::vector<int>& todecompress)
{
    std::vector<int> tdc = todecompress;
    
    int len = todecompress.size();
    int dlen = todecompress[len-1];

    todecompress = std::vector<int>(dlen, 0);

    int di = 0;
    for (int i = 0; i < len-1; i++)
    {
        if (tdc[i] > 0) // there are no 0 values in the compressed vector
        {
            todecompress[di] = tdc[i];
            di++;
        }
        else
            di += std::abs(tdc[i]);
    }
}

int gentools::getmaxdim(std::vector<std::vector<int>>* elementlist)
{
    int highestdim = -1;
    for (int i = 0; i < 8; i++)
    {
        element el(i);
        int eldim = el.getelementdimension();
        if (elementlist->at(i).size() > 0 && eldim > highestdim)
            highestdim = eldim;
    }
    return highestdim;
}

int gentools::sum(std::vector<int>& tosum)
{
    int out = 0;
    for (int i = 0; i < tosum.size(); i++)
        out += tosum[i];
    return out;
}

double gentools::sum(std::vector<double>& tosum)
{
    double out = 0;
    for (int i = 0; i < tosum.size(); i++)
        out += tosum[i];
    return out;
}

void gentools::splitatcolon(std::string tosplit, std::string& first, std::string& last)
{
    int colonpos = -1;
    for (int i = 0; i < tosplit.size(); i++)
    {
        if (tosplit[i] == ':')
        {
            colonpos = i;
            break;
        }
    }
    
    if (colonpos == -1)
    {
        first = "";
        last = tosplit;
    }
    else
    {
        first = tosplit.substr(0, colonpos);
        last = tosplit.substr(colonpos+1, tosplit.size()-colonpos-1);
    }
}

void gentools::findtruefalse(std::vector<bool>& invec, indexmat& trueinds, indexmat& falseinds, std::vector<int>& renum)
{
    int numtot = invec.size();
    int numtrue = counttrue(invec);

    renum = std::vector<int>(numtot);

    trueinds = indexmat(numtrue, 1);
    falseinds = indexmat(numtot-numtrue, 1);
    int* tiptr = trueinds.getvalues();
    int* fiptr = falseinds.getvalues();

    int ti = 0, fi = 0;
    for (int i = 0; i < numtot; i++)
    {
        if (invec[i])
        {
            renum[i] = ti;
            tiptr[ti] = i;
            ti++;
        }
        else
        {
            renum[i] = fi;
            fiptr[fi] = i;
            fi++;
        }
    }
}

void gentools::inoutorient(int physreg, std::vector<bool>& flipit)
{
    elements* els = universe::getrawmesh()->getelements();
    disjointregions* drs = universe::getrawmesh()->getdisjointregions();
    physicalregions* prs = universe::getrawmesh()->getphysicalregions();
    int totnumedges = els->count(1);
    
    std::vector<int> edgeinfo(totnumedges, 0);
    
    std::vector<int> ders = prs->get(physreg)->getdisjointregions(1);
    
    int numedgesinpr = 0;
    for (int i = 0; i < ders.size(); i++)
    {
        int rb = drs->getrangebegin(ders[i]);
        int ne = drs->countelements(ders[i]);
        
        for (int j = 0; j < ne; j++)
            edgeinfo[rb+j] = 2;
        
        numedgesinpr += ne;
    }
    
    // This loop on all edges is needed in case the physical region is made up of multiple disconnected islands:
    for (int i = 0; i < ders.size(); i++)
    {
        int rb = drs->getrangebegin(ders[i]);
        int ne = drs->countelements(ders[i]);
        
        for (int j = 0; j < ne; j++)
        {
            int firstnode = els->getsubelement(0, 1, rb+j, 0);
            // Process the whole island:
            inoutorient(firstnode, edgeinfo, true, false);
        }
    }
    
    flipit = std::vector<bool>(numedgesinpr);
        
    int index = 0;
    for (int i = 0; i < ders.size(); i++)
    {
        int rb = drs->getrangebegin(ders[i]);
        int ne = drs->countelements(ders[i]);
        
        for (int j = 0; j < ne; j++)
        {
            if (edgeinfo[rb+j] == 1)
                flipit[index] = false;
            else
                flipit[index] = true;
            index++;
        }
    }
}

// For edge number i value 'edgestatus[i]' is -1 if it is flipped, 1 if it is not flipped, 2 if
// it is in the physical region but not yet processed and 0 if it is not in the physical region.
void gentools::inoutorient(int startnode, std::vector<int>& edgestatus, bool isoutward, bool isrecursivecall)
{   
    elements* els = universe::getrawmesh()->getelements();
    
    std::vector<int> eon = els->getedgesonnode(startnode);

    for (int e = 0; e < eon.size(); e++)
    {
        int ce = eon[e];
        int firstnode = els->getsubelement(0, 1, ce, 0);
        int lastnode = els->getsubelement(0, 1, ce, 1);
            
        if (edgestatus[ce] == 2)
        {
            if (isoutward)
            {
                if (firstnode == startnode)
                    edgestatus[ce] = 1;
                else
                    edgestatus[ce] = -1;
            }
            else
            {
                if (firstnode == startnode)
                    edgestatus[ce] = -1;
                else
                    edgestatus[ce] = 1;
            }
            
            if (firstnode == startnode)
                inoutorient(lastnode, edgestatus, not(isoutward), true);
            else
                inoutorient(firstnode, edgestatus, not(isoutward), true);
        }
        else
        {
            if (isrecursivecall && edgestatus[ce] != 0)
            {
                bool isedgeoutward = (firstnode == startnode && edgestatus[ce] == 1 || lastnode == startnode && edgestatus[ce] == -1);

                if (isoutward != isedgeoutward)
                {
                    std::cout << "Error in 'gentools' namespace: reorienting the edges to have them all pointing either inwards or outwards at every node is impossible on the requested physical region for the mesh provided" << std::endl;
                    abort();
                }
            }
        }
    }
}

void gentools::fixatoverlap(std::vector<std::vector<int>>& cellvalues)
{
    std::shared_ptr<dtracker> dt = universe::getrawmesh()->getdtracker();
    int numneighbours = dt->countneighbours();
    std::vector<int> myneighbours = dt->getneighbours();

    if (dt->isdefined() == false || dt->isoverlap() == false || slmpi::count() == 1)
        return;

    physicalregions* prs = universe::getrawmesh()->getphysicalregions();

    std::vector<std::vector<int>> cellvalsforneighbours(numneighbours);
    std::vector<std::vector<int>> cellvalsfromneighbours(numneighbours);

    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];

        physicalregion* iopr = prs->get(dt->getinneroverlap(cn));
        physicalregion* oopr = prs->get(dt->getouteroverlap(cn));

        cellvalsforneighbours[n].resize(iopr->countelements());
        cellvalsfromneighbours[n].resize(oopr->countelements());

        std::vector<std::vector<int>>* ellist = iopr->getelementlist();

        int index = 0;
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < ellist->at(i).size(); j++)
            {
                cellvalsforneighbours[n][index] = cellvalues[i][ellist->at(i)[j]];
                index++;
            }
        }
    }

    slmpi::exchange(myneighbours, cellvalsforneighbours, cellvalsfromneighbours);

    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];

        std::vector<std::vector<int>>* ellist = prs->get(dt->getouteroverlap(cn))->getelementlist();

        int index = 0;
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < ellist->at(i).size(); j++)
            {
                // The owner of the inner overlap decides:
                cellvalues[i][ellist->at(i)[j]] = cellvalsfromneighbours[n][index];
                index++;
            }
        }
    }
}

void gentools::getedgesininnerinterfaces(std::vector<std::vector<int>>& iiedgelists, std::vector<std::vector<int>>& oiedgelistspreallocated)
{
    elements* els = universe::getrawmesh()->getelements();
    physicalregions* prs = universe::getrawmesh()->getphysicalregions();

    std::shared_ptr<dtracker> dt = universe::getrawmesh()->getdtracker();
    int numneighbours = dt->countneighbours();
    std::vector<int> myneighbours = dt->getneighbours();

    iiedgelists = std::vector<std::vector<int>>(numneighbours, std::vector<int>(0));
    oiedgelistspreallocated = std::vector<std::vector<int>>(numneighbours, std::vector<int>(0));

    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];

        std::vector<int> innerinterfaces = {}, outerinterfaces = {};
        if (dt->isoverlap())
        {
            innerinterfaces = dt->getinneroverlapinterface(cn);
            outerinterfaces = dt->getouteroverlapinterface(cn);
        }
        else
            innerinterfaces = dt->getnooverlapinterface(cn);

        std::vector< std::vector<std::vector<int>>* > innerellists(innerinterfaces.size()), outerellists(outerinterfaces.size());
        for (int i = 0; i < innerinterfaces.size(); i++)
            innerellists[i] = prs->get(innerinterfaces[i])->getelementlist();
        for (int i = 0; i < outerinterfaces.size(); i++)
            outerellists[i] = prs->get(outerinterfaces[i])->getelementlist();

        std::vector<bool> isinii, isinoi;

        int numinneredges = els->istypeinelementlists(1, innerellists, isinii, false);
        int numouteredges = numinneredges;
        if (dt->isoverlap())
            numouteredges = els->istypeinelementlists(1, outerellists, isinoi, false);

        iiedgelists[n] = std::vector<int>(2*numinneredges, -1);
        oiedgelistspreallocated[n] = std::vector<int>(2*numouteredges, -1);

        int index = 0;
        for (int e = 0; e < isinii.size(); e++)
        {
            if (isinii[e])
            {
                iiedgelists[n][2*index] = e;
                index++;
            }
        }
    }
}

std::vector<int> gentools::appendneighbourvalues(std::vector<double>& toappendto, std::vector<double>& toappend, int togroup)
{
    std::shared_ptr<dtracker> dt = universe::getrawmesh()->getdtracker();
    int numneighbours = dt->countneighbours();
    std::vector<int> neighbours = dt->getneighbours();
    
    int origlen = toappendto.size();
    
    std::vector<int> grouped = {}, appendlens = {};
    if (dt->isdefined() && slmpi::count() > 1)
    {
        std::vector<int> snds(2*numneighbours);
        for (int n = 0; n < numneighbours; n++)
        {
            snds[2*n+0] = togroup;
            snds[2*n+1] = toappend.size();
        }
        slmpi::exchange(neighbours, snds, appendlens);

        grouped = extract(appendlens, 2, 0);
    }
    
    toappendto.resize(origlen + toappend.size() + sum(appendlens));
    
    for (int i = 0; i < toappend.size(); i++)
        toappendto[origlen+i] = toappend[i];

    if (dt->isdefined() && slmpi::count() > 1)
    {
        std::vector<int> sendlens(numneighbours, toappend.size());
        std::vector<std::vector<double>> snds(numneighbours, toappend);
        std::vector<double*> sendbuffers(numneighbours), receivebuffers(numneighbours);

        int index = origlen + toappend.size();
        for (int n = 0; n < numneighbours; n++)
        {
            sendbuffers[n] = snds[n].data();
            receivebuffers[n] = toappendto.data() + index;

            index += appendlens[n];
        }
        slmpi::exchange(neighbours, sendlens, sendbuffers, appendlens, receivebuffers);
    }

    return grouped;
}

std::vector<int> gentools::getactiveelements(std::vector<int> disjregs, std::vector<bool>& isnodeactive)
{
    elements* els = universe::getrawmesh()->getelements();
    disjointregions* drs = universe::getrawmesh()->getdisjointregions();
    
    std::vector<int> nns = {1,2,3,4,4,8,6,5};

    int cnt = 0;
    for (int i = 0; i < disjregs.size(); i++)
        cnt += drs->countelements(disjregs[i]);
    
    std::vector<int> elemnums(cnt);
    
    int index = 0;
    for (int i = 0; i < disjregs.size(); i++)
    {
        int tn = drs->getelementtypenumber(disjregs[i]);
        int ne = drs->countelements(disjregs[i]);
        int rb = drs->getrangebegin(disjregs[i]);
        
        for (int j = 0; j < ne; j++)
        {
            for (int n = 0; n < nns[tn]; n++)
            {
                int nd = els->getsubelement(0, tn, rb+j, n);
                if (isnodeactive[nd])
                {
                    elemnums[index] = rb+j;
                    index++;
                    
                    break;
                }
            }
        }
    }
    
    elemnums.resize(index);
    
    return elemnums;
}

