#include "rawdisk.h"


rawdisk::rawdisk(int physreg, std::shared_ptr<rawshape> centerpoint, double radius, int nummeshpts)
{
    if (centerpoint->getdimension() != 0)
    {
        std::cout << "Error in 'rawdisk' object: expected a point shape for the disk center" << std::endl;
        abort();
    }

    if (nummeshpts%4 != 0 || nummeshpts <= 0)
    {
        std::cout << "Error in 'rawdisk' object: the structured disk meshing code requires a number of mesh nodes multiple of 4" << std::endl;
        abort();
    }

    if (radius == 0)
    {
        std::cout << "Error in 'rawdisk' object: disk radius cannot be zero" << std::endl;
        abort();
    }

    myphysicalregion = physreg;

    mynummeshpoints = nummeshpts;

    mycenterpoint = centerpoint;

    myradius = radius;

    // Sons will be defined when meshing

    mesh();
}


std::shared_ptr<rawshape> rawdisk::extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection)
{
    return std::shared_ptr<rawextrusion>(new rawextrusion(physreg, shared_from_this(), height, numlayers, extrudedirection));
}

std::shared_ptr<rawshape> rawdisk::duplicate(void)
{
    std::shared_ptr<rawdisk> out(new rawdisk);
    *out = *this;

    out->sons = geotools::duplicate(sons);
    out->mycenterpoint = mycenterpoint->duplicate();

    out->replicatelinks(shared_from_this());

    return out;    
}

void rawdisk::setphysicalregion(int physreg)
{
    myphysicalregion = physreg;
}

int rawdisk::getdimension(void)
{
    return 2;
}

std::string rawdisk::getname(void)
{
    return "disk";
}


std::vector<std::shared_ptr<rawshape>> rawdisk::getsons(void) 
{ 
    return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawdisk::getsubshapes(void)
{
    return geotools::concatenate({sons,{mycenterpoint}});
}

void rawdisk::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    mycenterpoint = subshapes[subshapes.size()-1];
    subshapes.pop_back();
    sons = subshapes;
}

int rawdisk::getphysicalregion(void) 
{ 
    return myphysicalregion; 
}

std::vector<double>* rawdisk::getcoords(void) 
{ 
    return &mycoords; 
}

std::vector<std::vector<int>>* rawdisk::getelems(void) 
{ 
    return &myelems; 
}

std::shared_ptr<rawshape> rawdisk::getpointer(void) 
{ 
    return shared_from_this(); 
}


void rawdisk::mesh(void)
{
    // Get the coordinates of the center point:
    std::vector<double> centercoords = *(mycenterpoint->getcoords());

    // Define the center point:
    std::shared_ptr<rawshape> pc = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0], centercoords[1], centercoords[2]}));

    // Define the points for the 4 outer arcs:
    std::shared_ptr<rawshape> p1 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0]+myradius, centercoords[1], centercoords[2]}));
    std::shared_ptr<rawshape> p2 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0], centercoords[1]+myradius, centercoords[2]}));
    std::shared_ptr<rawshape> p3 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0]-myradius, centercoords[1], centercoords[2]}));
    std::shared_ptr<rawshape> p4 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0], centercoords[1]-myradius, centercoords[2]}));

    // Define the points for the 4 corners of the inner square:
    double squareratio = 0.5;
    
    std::shared_ptr<rawshape> p5 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0]+squareratio*myradius, centercoords[1], centercoords[2]}));
    std::shared_ptr<rawshape> p6 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0], centercoords[1]+squareratio*myradius, centercoords[2]}));
    std::shared_ptr<rawshape> p7 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0]-squareratio*myradius, centercoords[1], centercoords[2]}));
    std::shared_ptr<rawshape> p8 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0], centercoords[1]-squareratio*myradius, centercoords[2]}));


    // Create the disk contour (4 arcs):
    std::shared_ptr<rawshape> arc1 = std::shared_ptr<rawarc>(new rawarc(-1, {p1,p2,pc}, mynummeshpoints/4+1));
    std::shared_ptr<rawshape> arc2 = std::shared_ptr<rawarc>(new rawarc(-1, {p2,p3,pc}, mynummeshpoints/4+1));
    std::shared_ptr<rawshape> arc3 = std::shared_ptr<rawarc>(new rawarc(-1, {p3,p4,pc}, mynummeshpoints/4+1));
    std::shared_ptr<rawshape> arc4 = std::shared_ptr<rawarc>(new rawarc(-1, {p4,p1,pc}, mynummeshpoints/4+1));

    // Create the remaining lines:
    std::shared_ptr<rawshape> l1 = std::shared_ptr<rawline>(new rawline(-1, {p1,p5}, mynummeshpoints/4+1));
    std::shared_ptr<rawshape> l2 = std::shared_ptr<rawline>(new rawline(-1, {p2,p6}, mynummeshpoints/4+1));
    std::shared_ptr<rawshape> l3 = std::shared_ptr<rawline>(new rawline(-1, {p3,p7}, mynummeshpoints/4+1));
    std::shared_ptr<rawshape> l4 = std::shared_ptr<rawline>(new rawline(-1, {p4,p8}, mynummeshpoints/4+1));

    std::shared_ptr<rawshape> l5 = std::shared_ptr<rawline>(new rawline(-1, {p5,p6}, mynummeshpoints/4+1));
    std::shared_ptr<rawshape> l6 = std::shared_ptr<rawline>(new rawline(-1, {p6,p7}, mynummeshpoints/4+1));
    std::shared_ptr<rawshape> l7 = std::shared_ptr<rawline>(new rawline(-1, {p7,p8}, mynummeshpoints/4+1));
    std::shared_ptr<rawshape> l8 = std::shared_ptr<rawline>(new rawline(-1, {p8,p5}, mynummeshpoints/4+1));

    // Create the 5 quadrangles that make the disk:
    std::shared_ptr<rawshape> q1 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, {arc1,l2,l5,l1}));
    std::shared_ptr<rawshape> q2 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, {arc2,l3,l6,l2}));
    std::shared_ptr<rawshape> q3 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, {arc3,l4,l7,l3}));
    std::shared_ptr<rawshape> q4 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, {arc4,l1,l8,l4}));
    std::shared_ptr<rawshape> q5 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, {l5,l6,l7,l8}));



    ///// Group the mesh of the 5 quadrangles:

    sons = {arc1, arc2, arc3, arc4};

    mycoords = geotools::appendcoords({q1,q2,q3,q4,q5});
    myelems = geotools::appendelems({q1,q2,q3,q4,q5});
}



