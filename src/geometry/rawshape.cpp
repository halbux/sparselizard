#include "rawshape.h"


void rawshape::deform(expression xdeform, expression ydeform, expression zdeform)
{
	std::vector<double>* mycoords = getcoords();

	int numnodes = mycoords->size()/3;

	std::vector<double> xcoords(numnodes);
	std::vector<double> ycoords(numnodes);
	std::vector<double> zcoords(numnodes);

	for (int i = 0; i < numnodes; i++)
	{
		xcoords[i] = mycoords->at(3*i+0);
		ycoords[i] = mycoords->at(3*i+1);
		zcoords[i] = mycoords->at(3*i+2);
	}

	std::vector<double> xdeformvec = xdeform.evaluate(xcoords,ycoords,zcoords);
	std::vector<double> ydeformvec = ydeform.evaluate(xcoords,ycoords,zcoords);
	std::vector<double> zdeformvec = zdeform.evaluate(xcoords,ycoords,zcoords);

	for (int i = 0; i < numnodes; i++)
	{
		mycoords->at(3*i+0) += xdeformvec[i];
		mycoords->at(3*i+1) += ydeformvec[i];
		mycoords->at(3*i+2) += zdeformvec[i];
	}

	
	// Also deform the sons:
	std::vector<std::shared_ptr<rawshape>> sons = getsons();
	for (int i = 0; i < sons.size(); i++)
		sons[i]->deform(xdeform, ydeform, zdeform);
}

void rawshape::shift(double shiftx, double shifty, double shiftz)
{
	std::vector<double>* mycoords = getcoords();

	int numnodes = mycoords->size()/3;

	for (int i = 0; i < numnodes; i++)
	{
		mycoords->at(3*i+0) += shiftx;
		mycoords->at(3*i+1) += shifty;
		mycoords->at(3*i+2) += shiftz;
	}


	// Also shift the sons:
	std::vector<std::shared_ptr<rawshape>> sons = getsons();
	for (int i = 0; i < sons.size(); i++)
		sons[i]->shift(shiftx, shifty, shiftz);
}

void rawshape::rotate(double alphax, double alphay, double alphaz)
{
	// Get the coordinates of the shape to rotate:
	std::vector<double>* mycoords = getcoords();

    int numberofnodes = mycoords->size()/3;

    // Convert input degrees to radians:
    double pi = 3.1415926535897932384;
    alphax = alphax*2*pi/360;
    alphay = alphay*2*pi/360;
    alphaz = alphaz*2*pi/360;
    
	// Define the rotation matrix R = Rx*Ry*Rz:
	double Rxx = std::cos(alphay)*std::cos(alphaz); 
	double Rxy = -std::cos(alphay)*std::sin(alphaz); 
	double Rxz = std::sin(alphay); 
	double Ryx = std::cos(alphax)*std::sin(alphaz) + std::cos(alphaz)*std::sin(alphax)*std::sin(alphay); 
	double Ryy = std::cos(alphax)*std::cos(alphaz) - std::sin(alphax)*std::sin(alphay)*std::sin(alphaz); 
	double Ryz = -std::cos(alphay)*std::sin(alphax); 
	double Rzx = std::sin(alphax)*std::sin(alphaz) - std::cos(alphax)*std::cos(alphaz)*std::sin(alphay); 
	double Rzy = std::cos(alphaz)*std::sin(alphax) + std::cos(alphax)*std::sin(alphay)*std::sin(alphaz); 
	double Rzz = std::cos(alphax)*std::cos(alphay); 

	// Compute R*[coordx; coordy; coordz]:
	for (int nodenumber = 0; nodenumber < numberofnodes; nodenumber++)
	{
		double xcoord = mycoords->at(3*nodenumber+0);
		double ycoord = mycoords->at(3*nodenumber+1);
		double zcoord = mycoords->at(3*nodenumber+2);
		
		mycoords->at(3*nodenumber+0) = Rxx * xcoord + Rxy * ycoord + Rxz * zcoord;
		mycoords->at(3*nodenumber+1) = Ryx * xcoord + Ryy * ycoord + Ryz * zcoord;
		mycoords->at(3*nodenumber+2) = Rzx * xcoord + Rzy * ycoord + Rzz * zcoord;
	}


	// Also rotate the sons:
	std::vector<std::shared_ptr<rawshape>> sons = getsons();
	for (int i = 0; i < sons.size(); i++)
		sons[i]->rotate(alphax, alphay, alphaz);
}


std::shared_ptr<rawshape> rawshape::extrude(int physreg, double height, int numlayers)
{
	std::cout << "Error in 'rawshape' object: 'extrude' has not been defined for this shape" << std::endl;
	abort(); 
}

std::shared_ptr<rawshape> rawshape::duplicate(void)
{
	std::cout << "Error in 'rawshape' object: 'duplicate' has not been defined for this shape" << std::endl;
	abort(); 
}

void rawshape::flip(void)
{
	std::cout << "Error in 'rawshape' object: 'flip' has not been defined for this shape" << std::endl;
	abort(); 
}

void rawshape::setphysicalregion(int physreg) 
{
	std::cout << "Error in 'rawshape' object: 'setphysicalregion' has not been defined for this shape" << std::endl;
	abort(); 
}

int rawshape::getdimension(void)
{
	std::cout << "Error in 'rawshape' object: 'getdimension' has not been defined for this shape" << std::endl;
	abort(); 
}

std::string rawshape::getname(void)
{
	std::cout << "Error in 'rawshape' object: 'getname' has not been defined for this shape" << std::endl;
	abort(); 
}

std::vector<std::shared_ptr<rawshape>> rawshape::getsons(void)
{
	std::cout << "Error in 'rawshape' object: 'getsons' has not been defined for this shape" << std::endl;
	abort(); 
}

std::vector<std::shared_ptr<rawshape>> rawshape::getsubshapes(void)
{
	std::cout << "Error in 'rawshape' object: 'getsubshapes' has not been defined for this shape" << std::endl;
	abort(); 
}

int rawshape::getphysicalregion(void)
{
	std::cout << "Error in 'rawshape' object: 'getphysicalregion' has not been defined for this shape" << std::endl;
	abort(); 
}

std::vector<double>* rawshape::getcoords(void)
{
	std::cout << "Error in 'rawshape' object: 'getcoords' has not been defined for this shape" << std::endl;
	abort(); 
}

std::vector<std::vector<int>>* rawshape::getelems(void)
{
	std::cout << "Error in 'rawshape' object: 'getelems' has not been defined for this shape" << std::endl;
	abort(); 
}

std::shared_ptr<rawshape> rawshape::getpointer(void)
{
	std::cout << "Error in 'rawshape' object: 'getpointer' has not been defined for this shape" << std::endl;
	abort();
}

void rawshape::mesh(void)
{
	std::cout << "Error in 'rawshape' object: 'mesh' has not been defined for this shape" << std::endl;
	abort();
}



