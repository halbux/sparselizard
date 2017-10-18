#include "nodes.h"


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

void nodes::shift(double xshift, double yshift, double zshift)
{
    int numberofnodes = count();
    
	for (int nodenumber = 0; nodenumber < numberofnodes; nodenumber++)
	{
		mycoordinates[3*nodenumber+0] += xshift;
		mycoordinates[3*nodenumber+1] += yshift;
		mycoordinates[3*nodenumber+2] += zshift;
	}
}

void nodes::rotate(double alphax, double alphay, double alphaz)
{
    int numberofnodes = count();

    // Convert input degrees to radians:
    double pi = 3.1415926535897932384;
    alphax = alphax*2*pi/360;
    alphay = alphay*2*pi/360;
    alphaz = alphaz*2*pi/360;
    
	// Define the rotation matrix R = Rx*Ry*Rz:
	double Rxx = cos(alphay)*cos(alphaz); 
	double Rxy = -cos(alphay)*sin(alphaz); 
	double Rxz = sin(alphay); 
	double Ryx = cos(alphax)*sin(alphaz) + cos(alphaz)*sin(alphax)*sin(alphay); 
	double Ryy = cos(alphax)*cos(alphaz) - sin(alphax)*sin(alphay)*sin(alphaz); 
	double Ryz = -cos(alphay)*sin(alphax); 
	double Rzx = sin(alphax)*sin(alphaz) - cos(alphax)*cos(alphaz)*sin(alphay); 
	double Rzy = cos(alphaz)*sin(alphax) + cos(alphax)*sin(alphay)*sin(alphaz); 
	double Rzz = cos(alphax)*cos(alphay); 

	// Compute R*[coordx; coory; coordz]:
	for (int nodenumber = 0; nodenumber < numberofnodes; nodenumber++)
	{
		double xcoord = mycoordinates[3*nodenumber+0];
		double ycoord = mycoordinates[3*nodenumber+1];
		double zcoord = mycoordinates[3*nodenumber+2];
		
		mycoordinates[3*nodenumber+0] = Rxx * xcoord + Rxy * ycoord + Rxz * zcoord;
		mycoordinates[3*nodenumber+1] = Ryx * xcoord + Ryy * ycoord + Ryz * zcoord;
		mycoordinates[3*nodenumber+2] = Rzx * xcoord + Rzy * ycoord + Rzz * zcoord;
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

