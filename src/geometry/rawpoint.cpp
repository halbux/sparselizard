#include "rawpoint.h"


rawpoint::rawpoint(int physreg, std::vector<double> coords)
{
    myphysicalregion = physreg;
    mycoords = coords;
    myelems[0] = {0};

    if (coords.size() == 3)
        return;
    else
    {
        std::cout << "Error in 'rawpoint' object: expected a vector with 3 coordinates" << std::endl;
        abort();
    }
}

rawpoint::rawpoint(int physreg, std::vector<double>& allcoords, std::vector<std::vector<int>>& allelems)
{
    myphysicalregion = physreg;

    mycoords = allcoords;
    myelems = allelems;
}


std::shared_ptr<rawshape> rawpoint::extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection)
{
    return std::shared_ptr<rawextrusion>(new rawextrusion(physreg, shared_from_this(), height, numlayers, extrudedirection));
}

std::shared_ptr<rawshape> rawpoint::duplicate(void)
{
    std::shared_ptr<rawpoint> out(new rawpoint);
    *out = *this;

    out->sons = geotools::duplicate(sons);

    out->replicatelinks(shared_from_this());

    return out;    
}

void rawpoint::setphysicalregion(int physreg)
{
    myphysicalregion = physreg;
}

int rawpoint::getdimension(void)
{
    return 0;
}

std::string rawpoint::getname(void)
{
    return "point";
}


std::vector<std::shared_ptr<rawshape>> rawpoint::getsons(void) 
{ 
    return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawpoint::getsubshapes(void)
{
    return sons;
}

void rawpoint::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    sons = subshapes;
}

int rawpoint::getphysicalregion(void) 
{ 
    return myphysicalregion; 
}

std::vector<double>* rawpoint::getcoords(void) 
{ 
    return &mycoords; 
}

std::vector<std::vector<int>>* rawpoint::getelems(void) 
{ 
    return &myelems; 
}

std::shared_ptr<rawshape> rawpoint::getpointer(void) 
{ 
    return shared_from_this(); 
}



