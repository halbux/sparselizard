#include "rawarc.h"


rawarc::rawarc(int physreg, std::vector<std::shared_ptr<rawshape>> inputpoints, int nummeshpoints)
{
	if (inputpoints.size() != 3 || inputpoints[0]->getdimension() != 0 || inputpoints[1]->getdimension() != 0 || inputpoints[2]->getdimension() != 0)
	{
		std::cout << "Error in 'rawarc' object: expected three points in the arc definition (start, end, center)" << std::endl;
		abort();
	}

	if (nummeshpoints < 2)
	{
		std::cout << "Error in 'rawarc' object: expected at least two mesh nodes in the arc" << std::endl;
		abort();
	}

	myphysicalregion = physreg;

	mynummeshpoints = nummeshpoints;

	sons = inputpoints;

	mesh();
}
 

std::shared_ptr<rawshape> rawarc::extrude(int physreg, double height, int numlayers)
{
	return std::shared_ptr<rawextrusion>(new rawextrusion(physreg, sons, {shared_from_this()}, height, numlayers));
}

std::shared_ptr<rawshape> rawarc::duplicate(void)
{
    std::shared_ptr<rawarc> out(new rawarc);
    *out = *this;

	out->sons = geotools::duplicate(sons);

    return out;	
}

void rawarc::flip(void)
{
	sons = {sons[1],sons[0],sons[2]};////!!! extra rotationdir in private otherwise problem with flip: rot dir changes

    mycoords = geotools::flipcoords(mycoords);
}

void rawarc::setphysicalregion(int physreg)
{
	myphysicalregion = physreg;
}

int rawarc::getdimension(void)
{
	return 1;
}

std::string rawarc::getname(void)
{
	return "arc";
}


std::vector<std::shared_ptr<rawshape>> rawarc::getsons(void) 
{ 
	return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawarc::getsubshapes(void)
{
	return sons;
}

int rawarc::getphysicalregion(void) 
{ 
	return myphysicalregion; 
}

std::vector<double>* rawarc::getcoords(void) 
{ 
	return &mycoords; 
}

std::vector<std::vector<int>>* rawarc::getelems(void) 
{ 
	return &myelems; 
}

std::shared_ptr<rawshape> rawarc::getpointer(void) 
{ 
	return shared_from_this(); 
}


void rawarc::mesh(void)
{
	int numelems = mynummeshpoints-1;

	mycoords.resize(3*mynummeshpoints);
	myelems[1].resize(2*numelems);

	// p1 is the first point, p2 the last one and p3 the center point:
	std::vector<double>* p1coords = sons[0]->getcoords();
	std::vector<double>* p2coords = sons[1]->getcoords();
	std::vector<double>* p3coords = sons[2]->getcoords();

	// Give an error if not in a circle:
	double radius1 = std::sqrt(std::pow(p1coords->at(0)-p3coords->at(0),2.0)+std::pow(p1coords->at(1)-p3coords->at(1),2.0));
	double radius2 = std::sqrt(std::pow(p2coords->at(0)-p3coords->at(0),2.0)+std::pow(p2coords->at(1)-p3coords->at(1),2.0));

	if (radius1 == 0 || std::abs((radius1-radius2)/radius1) > 1e-10)
	{
		std::cout << "Error in 'rawarc' object: rawpoints 1 and 2 provided for arc should be on a circle with center rawpoint 3" << std::endl;
		abort();
	}
	double radius = radius1;

	// Get the angle with the x axis of the first and last point on the arc:
	double angle1 = std::acos((p1coords->at(0)-p3coords->at(0))/radius);
	double angle2 = std::acos((p2coords->at(0)-p3coords->at(0))/radius);
	if (p1coords->at(1) < p3coords->at(1))
		angle1 = -angle1;
	if (p2coords->at(1) < p3coords->at(1))
		angle2 = -angle2;
	// Angle between two consecutive points in the mesh:
	double deltaangle = (angle2-angle1)/numelems;

	for (int i = 0; i < mynummeshpoints; i++)
	{
		mycoords[3*i+0] = p3coords->at(0) + radius*cos(angle1+i*deltaangle);
		mycoords[3*i+1] = p3coords->at(1) + radius*sin(angle1+i*deltaangle);;
		mycoords[3*i+2] = p3coords->at(2);
	}


	for (int i = 0; i < numelems; i++)
	{
		myelems[1][2*i+0] = i;
		myelems[1][2*i+1] = (i+1)%mynummeshpoints;
	}
}

