#include "nodes.h"
#include "universe.h"
#include "gentools.h"


nodes::nodes(void) {}

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

std::vector<int> nodes::removeduplicates(void)
{
    // 'noderenumbering' will give the renumbering corresponding to removed duplicates:
    std::vector<int> noderenumbering;
    int numberofnonduplicates = gentools::removeduplicates(mycoordinates, noderenumbering);

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

std::vector<double> nodes::getgeometrydimension(void)
{
    std::vector<double> bnds = gentools::getcoordbounds(mycoordinates);

    return {std::abs(bnds[0]-bnds[1]), std::abs(bnds[2]-bnds[3]), std::abs(bnds[4]-bnds[5])};
}

std::vector<double> nodes::getnoisethreshold(void)
{
    double rons = universe::roundoffnoiselevel;
    
    std::vector<double> gds = getgeometrydimension();
    
    return std::vector<double> {rons*gds[0], rons*gds[1], rons*gds[2]};
}

void nodes::fixifaxisymmetric(void)
{
    if (universe::isaxisymmetric == false)
        return;

    double xnoiselevel = getnoisethreshold()[0];

    int numberofnodes = count();
    for (int i = 0; i < numberofnodes; i++)
    {
        double curx = mycoordinates[3*i+0];
        if (curx < 0)
        {
            if (std::abs(curx) < xnoiselevel)
                mycoordinates[3*i+0] = 0.0;
            else
            {
                logs log;
                log.msg() << "Error in 'nodes' object: expected only positive x node coordinates with axisymmetry (found a node at x = " << curx << " which is out of noise range)" << std::endl;
                log.error();
            }
        }
    }
}

