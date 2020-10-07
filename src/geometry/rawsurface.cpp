#include "rawsurface.h"


rawsurface::rawsurface(int physreg, std::vector<double>& allcoords, std::vector<std::vector<int>>& allelems)
{
    myphysicalregion = physreg;

    mycoords = allcoords;
    myelems = allelems;
}



std::shared_ptr<rawshape> rawsurface::extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection)
{
    return std::shared_ptr<rawextrusion>(new rawextrusion(physreg, shared_from_this(), height, numlayers, extrudedirection));
}

std::shared_ptr<rawshape> rawsurface::duplicate(void)
{
    std::shared_ptr<rawsurface> out(new rawsurface);
    *out = *this;

    out->sons = geotools::duplicate(sons);

    out->replicatelinks(shared_from_this());

    return out;    
}

void rawsurface::setphysicalregion(int physreg)
{
    myphysicalregion = physreg;
}

int rawsurface::getdimension(void)
{
    return 2;
}

std::string rawsurface::getname(void)
{
    return "surface";
}


std::vector<std::shared_ptr<rawshape>> rawsurface::getsons(void) 
{ 
    return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawsurface::getsubshapes(void)
{
    return sons;
}

void rawsurface::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    sons = subshapes;
}

int rawsurface::getphysicalregion(void) 
{ 
    return myphysicalregion; 
}

std::vector<double>* rawsurface::getcoords(void) 
{ 
    return &mycoords; 
}

std::vector<std::vector<int>>* rawsurface::getelems(void) 
{ 
    return &myelems; 
}

std::shared_ptr<rawshape> rawsurface::getpointer(void) 
{ 
    return shared_from_this(); 
}


void rawsurface::mesh(void)
{

}



