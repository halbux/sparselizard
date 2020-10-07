#include "rawunion.h"


rawunion::rawunion(int physreg, std::vector<std::shared_ptr<rawshape>> input)
{
    sons = input;
    
    if (input.size() == 0)
    {
        std::cout << "Error in 'rawunion' object: expected at least one shape to unite" << std::endl;
        abort();
    }
    
    int mydimension = input[0]->getdimension();
    for (int i = 1; i < input.size(); i++)
    {
        if (input[i]->getdimension() != mydimension)
        {
            std::cout << "Error in 'rawunion' object: all shapes must have the same dimension" << std::endl;
            abort();
        }
    }

    myphysicalregion = physreg;
        
    mesh();
}

std::shared_ptr<rawshape> rawunion::extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection)
{
    std::vector<std::shared_ptr<rawshape>> extrudedshapes(sons.size());
    for (int i = 0; i < sons.size(); i++)
        extrudedshapes[i] = sons[i]->extrude(physreg, height, numlayers, extrudedirection);

    return std::shared_ptr<rawunion>(new rawunion(physreg, extrudedshapes));
}

std::shared_ptr<rawshape> rawunion::duplicate(void)
{
    std::shared_ptr<rawunion> out(new rawunion);
    *out = *this;

    out->sons = geotools::duplicate(sons);

    out->replicatelinks(shared_from_this());

    return out;    
}

void rawunion::setphysicalregion(int physreg)
{
    myphysicalregion = physreg;
}

int rawunion::getdimension(void)
{
    return sons[0]->getdimension();
}

std::string rawunion::getname(void)
{
    return "union";
}


std::vector<std::shared_ptr<rawshape>> rawunion::getsons(void) 
{ 
    return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawunion::getsubshapes(void)
{
    return sons;
}

void rawunion::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    sons = subshapes;
}

int rawunion::getphysicalregion(void) 
{ 
    return myphysicalregion; 
}

std::vector<double>* rawunion::getcoords(void) 
{ 
    return &mycoords; 
}

std::vector<std::vector<int>>* rawunion::getelems(void) 
{ 
    return &myelems; 
}

std::shared_ptr<rawshape> rawunion::getpointer(void) 
{ 
    return shared_from_this(); 
}


void rawunion::mesh(void)
{
    // Append the coordinate of the building blocks:
    mycoords = geotools::appendcoords(sons);
    // Append the elements of the building blocks:
    myelems = geotools::appendelems(sons);
}



