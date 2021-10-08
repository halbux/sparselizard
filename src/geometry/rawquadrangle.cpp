#include "rawquadrangle.h"


rawquadrangle::rawquadrangle(int physreg, std::vector<std::shared_ptr<rawshape>> inputpoints, std::vector<int> nummeshpoints)
{
    if (inputpoints.size() != 4 || inputpoints[0]->getdimension() != 0 || inputpoints[1]->getdimension() != 0 || inputpoints[2]->getdimension() != 0 || inputpoints[3]->getdimension() != 0)
    {
        std::cout << "Error in 'rawquadrangle' object: expected four points in the quadrangle definition" << std::endl;
        abort();
    }

    if (nummeshpoints.size() != 4)
    {
        std::cout << "Error in 'rawquadrangle' object: expected a vector of length four to define the number of mesh nodes" << std::endl;
        abort();
    }

    std::vector<std::shared_ptr<rawshape>> lns(4);
    
    for (int i = 0; i < 4; i++)
        lns[i] = std::shared_ptr<rawline>(new rawline(-1, {inputpoints[i], inputpoints[(i+1)%4]}, nummeshpoints[i]));


    myphysicalregion = physreg;

    sons = lns;

    mesh();
}

rawquadrangle::rawquadrangle(int physreg, std::vector<std::shared_ptr<rawshape>> inputlines)
{
    if (inputlines.size() != 4 || inputlines[0]->getdimension() != 1 || inputlines[1]->getdimension() != 1 || inputlines[2]->getdimension() != 1 || inputlines[3]->getdimension() != 1)
    {
        std::cout << "Error in 'rawquadrangle' object: expected four lines in the quadrangle definition" << std::endl;
        abort();
    }

    myphysicalregion = physreg;

    sons = geotools::orient(inputlines);

    mesh();
}



std::shared_ptr<rawshape> rawquadrangle::extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection)
{
    return std::shared_ptr<rawextrusion>(new rawextrusion(physreg, shared_from_this(), height, numlayers, extrudedirection));
}

std::shared_ptr<rawshape> rawquadrangle::duplicate(void)
{
    std::shared_ptr<rawquadrangle> out(new rawquadrangle);
    *out = *this;

    out->sons = geotools::duplicate(sons);

    out->replicatelinks(shared_from_this());

    return out;    
}

void rawquadrangle::setphysicalregion(int physreg)
{
    myphysicalregion = physreg;
}

int rawquadrangle::getdimension(void)
{
    return 2;
}

std::string rawquadrangle::getname(void)
{
    return "quadrangle";
}


std::vector<std::shared_ptr<rawshape>> rawquadrangle::getsons(void) 
{ 
    return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawquadrangle::getsubshapes(void)
{
    return sons;
}

void rawquadrangle::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    sons = subshapes;
}

int rawquadrangle::getphysicalregion(void) 
{ 
    return myphysicalregion; 
}

std::vector<double>* rawquadrangle::getcoords(void) 
{ 
    return &mycoords; 
}

std::vector<std::vector<int>>* rawquadrangle::getelems(void) 
{ 
    return &myelems; 
}

std::shared_ptr<rawshape> rawquadrangle::getpointer(void) 
{ 
    return shared_from_this(); 
}


void rawquadrangle::mesh(void)
{
    // Get the node coordinates on the lines:
    std::vector<std::vector<double>*> linescoords(4);
    std::vector<int> nummeshpts(4);

    for (int i = 0; i < 4; i++)
    {
        linescoords[i] = sons[i]->getcoords();
        nummeshpts[i] = linescoords[i]->size()/3;
    }

    // Give an error if lines facing each other do not have the same number of mesh nodes:
    if (nummeshpts[0] != nummeshpts[2] || nummeshpts[1] != nummeshpts[3])
    {
        std::cout << "Error in 'rawquadrangle' object: the number of nodes on edges facing each other should be equal" << std::endl;
        abort(); 
    }

    // Preallocate the coords and elems containers:
    mycoords.resize(3*nummeshpts[0]*nummeshpts[1]);
    myelems[3].resize(4*(nummeshpts[0]-1)*(nummeshpts[1]-1));

    // Get the coordinates of the corner nodes:
    std::vector<double> xcorner = {linescoords[0]->at(0), linescoords[1]->at(0), linescoords[2]->at(0), linescoords[3]->at(0)};
    std::vector<double> ycorner = {linescoords[0]->at(1), linescoords[1]->at(1), linescoords[2]->at(1), linescoords[3]->at(1)};
    std::vector<double> zcorner = {linescoords[0]->at(2), linescoords[1]->at(2), linescoords[2]->at(2), linescoords[3]->at(2)};

    
    int nki = nummeshpts[0];
    int neta = nummeshpts[1];

    // Loop on all layers in the ki direction:
    for (int i = 0; i < nki; i++)
    {    
        double mu = (double)i/(nki-1);

        // Coordinates of the first and last nodes in the line linking 
        // the current node in line 0 and its counterpart in line 2:
        double x1 = linescoords[0]->at(3*i+0);
        double xn = linescoords[2]->at(3*(nki-1-i)+0);
        double y1 = linescoords[0]->at(3*i+1);
        double yn = linescoords[2]->at(3*(nki-1-i)+1);
        double z1 = linescoords[0]->at(3*i+2);
        double zn = linescoords[2]->at(3*(nki-1-i)+2);

        for (int j = 0; j < neta; j++)
        {
            double lambda = (double)j/(neta-1);

            // Coordinates of the current node in the undeformed quadrangle:
            double xundef = (1-lambda)*((1-mu)*xcorner[0]+mu*xcorner[1]) + lambda*((1-mu)*xcorner[3]+mu*xcorner[2]);
            double yundef = (1-lambda)*((1-mu)*ycorner[0]+mu*ycorner[1]) + lambda*((1-mu)*ycorner[3]+mu*ycorner[2]);
            double zundef = (1-lambda)*((1-mu)*zcorner[0]+mu*zcorner[1]) + lambda*((1-mu)*zcorner[3]+mu*zcorner[2]);

            mycoords[3*i*neta+3*j + 0] = (1-lambda)*x1 + lambda*xn  +  mu*linescoords[1]->at(3*j+0) + (1-mu)*linescoords[3]->at(3*(neta-1-j)+0) - xundef;
            mycoords[3*i*neta+3*j + 1] = (1-lambda)*y1 + lambda*yn  +  mu*linescoords[1]->at(3*j+1) + (1-mu)*linescoords[3]->at(3*(neta-1-j)+1) - yundef;
            mycoords[3*i*neta+3*j + 2] = (1-lambda)*z1 + lambda*zn  +  mu*linescoords[1]->at(3*j+2) + (1-mu)*linescoords[3]->at(3*(neta-1-j)+2) - zundef;
        }
    }


    int currentelem = 0;
    for (int i = 0; i < nki-1; i++)
    {
        for (int j = 0; j < neta-1; j++)
        {
            myelems[3][4*currentelem+0] = i*neta+j;
            myelems[3][4*currentelem+1] = (i+1)*neta+j;
            myelems[3][4*currentelem+2] = (i+1)*neta+(j+1);
            myelems[3][4*currentelem+3] = i*neta+j+1;

            currentelem++;
        }
    }
}



