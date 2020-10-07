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

    sons = {inputpoints[0], inputpoints[1]};
    mycenterpoint = inputpoints[2];

    mesh();
}
 

std::shared_ptr<rawshape> rawarc::extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection)
{
    return std::shared_ptr<rawextrusion>(new rawextrusion(physreg, shared_from_this(), height, numlayers, extrudedirection));
}

std::shared_ptr<rawshape> rawarc::duplicate(void)
{
    std::shared_ptr<rawarc> out(new rawarc);
    *out = *this;

    out->sons = geotools::duplicate(sons);
    out->mycenterpoint = geotools::duplicate({mycenterpoint})[0];

    out->replicatelinks(shared_from_this());

    return out;    
}

void rawarc::flip(void)
{
    sons = {sons[1],sons[0]};

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
    return geotools::concatenate({sons,{mycenterpoint}});
}

void rawarc::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    mycenterpoint = subshapes[subshapes.size()-1];
    subshapes.pop_back();
    sons = subshapes;
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
    double pi = 3.1415926535897932384;
    
    int numelems = mynummeshpoints-1;

    mycoords.resize(3*mynummeshpoints);
    myelems[1].resize(2*numelems);

    // p1 is the first point, p2 the last one and pc the center point:
    std::vector<double> p1coords = *(sons[0]->getcoords());
    std::vector<double> p2coords = *(sons[1]->getcoords());
    std::vector<double> pccoords = *(mycenterpoint->getcoords());

    // Give an error if not in a circle:
    double radius1 = geotools::getdistance(p1coords, pccoords);
    double radius2 = geotools::getdistance(p2coords, pccoords);

    if (radius1 == 0.0 || std::abs((radius1-radius2)/radius1) > 1e-10)
    {
        std::cout << "Error in 'rawarc' object: rawpoints 1 and 2 provided for arc should be on a circle with center rawpoint 3" << std::endl;
        abort();
    }
    double radius = radius1;


    ///// Rotate the arc plane to have it parallel to the xy plane:
    double xaxisrot = geotools::getplanerotation("xrot", pccoords, p1coords, p2coords);
    if (xaxisrot != 0)
    {
        sons[0]->rotate(xaxisrot,0,0);
        sons[1]->rotate(xaxisrot,0,0);
        mycenterpoint->rotate(xaxisrot,0,0);

        p1coords = *(sons[0]->getcoords());
        p2coords = *(sons[1]->getcoords());
        pccoords = *(mycenterpoint->getcoords());
    }
    double yaxisrot = geotools::getplanerotation("yrot", pccoords, p1coords, p2coords);
    if (yaxisrot != 0)
    {
        sons[0]->rotate(0,yaxisrot,0);
        sons[1]->rotate(0,yaxisrot,0);
        mycenterpoint->rotate(0,yaxisrot,0);

        p1coords = *(sons[0]->getcoords());
        p2coords = *(sons[1]->getcoords());
        pccoords = *(mycenterpoint->getcoords());
    }
    ///// Rotated


    // Get the angle with the x axis of the first and last point on the arc:
    double angle1 = geotools::acos((p1coords[0]-pccoords[0])/radius);
    double angle2 = geotools::acos((p2coords[0]-pccoords[0])/radius);
    if (p1coords[1] < pccoords[1])
        angle1 = -angle1;
    if (p2coords[1] < pccoords[1])
        angle2 = -angle2;
    // Angle between two consecutive points in the mesh:
    double deltaangle = (angle2-angle1)/numelems;
    if (angle2-angle1 < 1e-10)
        deltaangle = (2*pi+angle2-angle1)/numelems;

    for (int i = 0; i < mynummeshpoints; i++)
    {
        mycoords[3*i+0] = pccoords[0] + radius*cos(angle1+i*deltaangle);
        mycoords[3*i+1] = pccoords[1] + radius*sin(angle1+i*deltaangle);
        mycoords[3*i+2] = pccoords[2];
    }


    for (int i = 0; i < numelems; i++)
    {
        myelems[1][2*i+0] = i;
        myelems[1][2*i+1] = i+1;
    }


    ///// Rotate everything back:
    if (xaxisrot != 0 || yaxisrot != 0)
    {
        rotate(0,-yaxisrot,0);
        rotate(-xaxisrot,0,0);
    }
}

