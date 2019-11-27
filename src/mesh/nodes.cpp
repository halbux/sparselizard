#include "nodes.h"
#include "geotools.h"


nodes::nodes(void) {}

nodes::nodes(double roundoffnoise)
{
    roundoffnoiselevel = roundoffnoise;
}

void nodes::setnumber(int numberofnodes)
{
    mycoordinates.resize(3*numberofnodes);
}

int nodes::count(void)
{
    return mycoordinates.size()/3;
}

std::vector<double>* nodes::getcoordinates(void)
{
    return &mycoordinates;
}

void nodes::shift(int physreg, double xshift, double yshift, double zshift)
{
    int numberofnodes = count();
    
    std::vector<bool> toshift;

    if (physreg < 0)
        toshift = std::vector<bool>(numberofnodes, true);
    else
    {
        toshift = std::vector<bool>(numberofnodes, false);
    
        disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();
        elements* myelements = universe::mymesh->getelements();
    
        // Get only the disjoint regions with highest dimension elements:
        std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions();
    
        for (int i = 0; i < selecteddisjregs.size(); i++)
        {
            int disjreg = selecteddisjregs[i];
            int numelems = mydisjointregions->countelements(disjreg);
            int elemtypenum = mydisjointregions->getelementtypenumber(disjreg);
            int rangebegin = mydisjointregions->getrangebegin(disjreg);
            int curvatureorder = myelements->getcurvatureorder();
            
            element myelem(elemtypenum, curvatureorder);

            for (int e = 0; e < numelems; e++)
            {
                for (int n = 0; n < myelem.countcurvednodes(); n++)
                    toshift[myelements->getsubelement(0, elemtypenum, rangebegin+e, n)] = true;
            }
        }
    }
    
    for (int n = 0; n < numberofnodes; n++)
    {
        if (toshift[n])
        {
            mycoordinates[3*n+0] += xshift;
            mycoordinates[3*n+1] += yshift;
            mycoordinates[3*n+2] += zshift;
        }
    }
}

void nodes::rotate(int physreg, double alphax, double alphay, double alphaz)
{
    int numberofnodes = count();
    
    std::vector<bool> torotate;

    if (physreg < 0)
        torotate = std::vector<bool>(numberofnodes, true);
    else
    {
        torotate = std::vector<bool>(numberofnodes, false);
    
        disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();
        elements* myelements = universe::mymesh->getelements();
    
        // Get only the disjoint regions with highest dimension elements:
        std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions();
    
        for (int i = 0; i < selecteddisjregs.size(); i++)
        {
            int disjreg = selecteddisjregs[i];
            int numelems = mydisjointregions->countelements(disjreg);
            int elemtypenum = mydisjointregions->getelementtypenumber(disjreg);
            int rangebegin = mydisjointregions->getrangebegin(disjreg);
            int curvatureorder = myelements->getcurvatureorder();
            
            element myelem(elemtypenum, curvatureorder);

            for (int e = 0; e < numelems; e++)
            {
                for (int n = 0; n < myelem.countcurvednodes(); n++)
                    torotate[myelements->getsubelement(0, elemtypenum, rangebegin+e, n)] = true;
            }
        }
    }
    
    std::vector<double> rotated = mycoordinates;
    geotools::rotate(alphax, alphay, alphaz, &rotated);
    
    for (int n = 0; n < numberofnodes; n++)
    {
        if (torotate[n])
        {
            mycoordinates[3*n+0] = rotated[3*n+0];
            mycoordinates[3*n+1] = rotated[3*n+1];
            mycoordinates[3*n+2] = rotated[3*n+2];
        }
    }
}

void nodes::scale(int physreg, double xscale, double yscale, double zscale)
{
    int numberofnodes = count();
    
    std::vector<bool> toscale;

    if (physreg < 0)
        toscale = std::vector<bool>(numberofnodes, true);
    else
    {
        toscale = std::vector<bool>(numberofnodes, false);
    
        disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();
        elements* myelements = universe::mymesh->getelements();
    
        // Get only the disjoint regions with highest dimension elements:
        std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions();
    
        for (int i = 0; i < selecteddisjregs.size(); i++)
        {
            int disjreg = selecteddisjregs[i];
            int numelems = mydisjointregions->countelements(disjreg);
            int elemtypenum = mydisjointregions->getelementtypenumber(disjreg);
            int rangebegin = mydisjointregions->getrangebegin(disjreg);
            int curvatureorder = myelements->getcurvatureorder();
            
            element myelem(elemtypenum, curvatureorder);

            for (int e = 0; e < numelems; e++)
            {
                for (int n = 0; n < myelem.countcurvednodes(); n++)
                    toscale[myelements->getsubelement(0, elemtypenum, rangebegin+e, n)] = true;
                    
            }
        }
    }
    
    for (int n = 0; n < numberofnodes; n++)
    {
        if (toscale[n])
        {
            mycoordinates[3*n+0] *= xscale;
            mycoordinates[3*n+1] *= yscale;
            mycoordinates[3*n+2] *= zscale;
        }
    }
}

void nodes::print(void)
{
    std::cout << "Number of nodes: " << count() << std::endl;
    std::cout << std::endl << "node | x coord | y coord | z coord" << std::endl << std::endl;
    
    int oldprecision = std::cout.precision();
    std::cout.precision(17);
    
    for (int i = 0; i < count(); i++)
        std::cout << std::setw(10) << std::left << i << std::setw(26) << std::left << mycoordinates[3*i+0] << std::setw(26) << std::left << mycoordinates[3*i+1] << std::setw(26) << std::left << mycoordinates[3*i+2] << std::endl;
    std::cout << std::endl;
    
    std::cout.precision(oldprecision);
}

std::vector<int> nodes::sortbycoordinates(void)
{
    int numberofnodes = count();
    
    // 'reorderingvector' gives the relation between the indexes before and after node sorting:
    std::vector<int> reorderingvector;
    myalgorithm::stablecoordinatesort(getnoisethreshold(), mycoordinates, reorderingvector);

    // Reorder the nodes.:
    reorder(reorderingvector);
    
    // sortedcoordinates = coordinates(reorderingvector,:).
    // sortedcoordinates(renumberingvector,:) = coordinates.
    std::vector<int> renumberingvector(numberofnodes);
    for (int i = 0; i < numberofnodes; i++)
        renumberingvector[reorderingvector[i]] = i;
    
    return renumberingvector;
}

// 'removeduplicates' assumes that the nodes are already sorted coordinatewise 
// thus the duplicated nodes follow each other in 'mycoordinates'.
std::vector<int> nodes::removeduplicates(void)
{
    int numberofnodes = count();
    
    // 'noderenumbering' will give the renumbering corresponding to removed duplicates:
    std::vector<int> noderenumbering;
    int numberofnonduplicates = myalgorithm::removeduplicatedcoordinates(getnoisethreshold(), mycoordinates, noderenumbering);

    for (int i = 0; i < noderenumbering.size(); i++)
    {
        if (noderenumbering[i] != i)
        {
            mycoordinates[3*noderenumbering[i]+0] = mycoordinates[3*i+0];
            mycoordinates[3*noderenumbering[i]+1] = mycoordinates[3*i+1];
            mycoordinates[3*noderenumbering[i]+2] = mycoordinates[3*i+2];
        }
    }
    // Remove unused space:
    mycoordinates.resize(numberofnonduplicates*3);
    
    return noderenumbering;
}

void nodes::reorder(std::vector<int>& nodereordering)
{
    int numberofnodes = count();
    
    // Update 'mycoordinates':
    std::vector<double> nodecoordinatescopy = mycoordinates;
    for (int i = 0; i < numberofnodes; i++)
    {
        mycoordinates[3*i+0] = nodecoordinatescopy[3*nodereordering[i]+0]; 
        mycoordinates[3*i+1] = nodecoordinatescopy[3*nodereordering[i]+1]; 
        mycoordinates[3*i+2] = nodecoordinatescopy[3*nodereordering[i]+2]; 
    }
}

double nodes::getgeometrydimension(int coord)
{
    int numberofnodes = count();
    
    double maxcoord = mycoordinates[3*0+coord], mincoord = mycoordinates[3*0+coord];
    for (int i = 1; i < numberofnodes; i++)
    {
        if (mycoordinates[3*i+coord] < mincoord)
            mincoord = mycoordinates[3*i+coord];
        if (mycoordinates[3*i+coord] > maxcoord)
            maxcoord = mycoordinates[3*i+coord];
    }
    return std::abs(maxcoord-mincoord);
}

std::vector<double> nodes::getnoisethreshold(void)
{
    return std::vector<double> {roundoffnoiselevel*getgeometrydimension(0), roundoffnoiselevel*getgeometrydimension(1), roundoffnoiselevel*getgeometrydimension(2)};
}

