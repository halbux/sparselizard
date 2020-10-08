#include "shape.h"

#include "rawpoint.h"
#include "rawline.h"
#include "rawarc.h"
#include "rawquadrangle.h"
#include "rawtriangle.h"
#include "rawdisk.h"
#include "rawsurface.h"
#include "rawextrusion.h"
#include "rawvolume.h"
#include "rawunion.h"


void shape::errornullpointer(void)
{
    if (rawshapeptr != NULL)
        return;
    else
    {
        std::cout << "Error in 'shape' object: cannot perform operation (shape is undefined)" << std::endl;
        abort();
    }
}

shape::shape(void) {}

shape::shape(std::shared_ptr<rawshape> inputptr)
{
    rawshapeptr = inputptr;
}

shape::shape(std::string shapename, int physreg, std::vector<double> coords)
{
    if (coords.size()%3 != 0)
    {
        std::cout << "Error in 'shape' object: length of coordinate vector for shape " << shapename << " should be a multiple of three" << std::endl;
        abort();
    }

    if (shapename == "point")
    {
        rawshapeptr = std::shared_ptr<rawpoint>(new rawpoint(physreg, coords));
        return;
    }
    if (shapename == "line")
    {
        rawshapeptr = std::shared_ptr<rawline>(new rawline(physreg, coords));
        return;
    }

    std::cout << "Error in 'shape' object: shape " << shapename << " does not accept this constructor or is unknown (try lower case)" << std::endl;
    abort();
}

shape::shape(std::string shapename, int physreg, std::vector<double> coords, int nummeshpts)
{
    std::vector<std::shared_ptr<rawshape>> cornerpts = geotools::coordstopoints(coords);

    if (shapename == "line")
    {
        rawshapeptr = std::shared_ptr<rawline>(new rawline(physreg, cornerpts, nummeshpts));
        return;
    }
    if (shapename == "arc")
    {
        rawshapeptr = std::shared_ptr<rawarc>(new rawarc(physreg, cornerpts, nummeshpts));
        return;
    }

    std::cout << "Error in 'shape' object: shape " << shapename << " does not accept this constructor or is unknown (try lower case)" << std::endl;
    abort();
}

shape::shape(std::string shapename, int physreg, std::vector<double> coords, std::vector<int> nummeshpts)
{
    std::vector<std::shared_ptr<rawshape>> cornerpts = geotools::coordstopoints(coords);

    if (shapename == "triangle")
    {
        rawshapeptr = std::shared_ptr<rawtriangle>(new rawtriangle(physreg, cornerpts, nummeshpts));
        return;
    }
    if (shapename == "quadrangle")
    {
        rawshapeptr = std::shared_ptr<rawquadrangle>(new rawquadrangle(physreg, cornerpts, nummeshpts));
        return;
    }

    std::cout << "Error in 'shape' object: shape " << shapename << " does not accept this constructor or is unknown (try lower case)" << std::endl;
    abort();
}

shape::shape(std::string shapename, int physreg, std::vector<shape> subshapes, int nummeshpts)
{
    // Get the rawshape pointer from all shapes:
    std::vector<std::shared_ptr<rawshape>> subrawshapes = geotools::getrawshapes(subshapes);

    if (shapename == "line")
    {
        rawshapeptr = std::shared_ptr<rawline>(new rawline(physreg, subrawshapes, nummeshpts));
        return;
    }
    if (shapename == "arc")
    {
        rawshapeptr = std::shared_ptr<rawarc>(new rawarc(physreg, subrawshapes, nummeshpts));
        return;
    }

    std::cout << "Error in 'shape' object: shape " << shapename << " does not accept this constructor or is unknown (try lower case)" << std::endl;
    abort();
}

shape::shape(std::string shapename, int physreg, std::vector<shape> subshapes, std::vector<int> nummeshpts)
{
    // Get the rawshape pointer from all shapes:
    std::vector<std::shared_ptr<rawshape>> subrawshapes = geotools::getrawshapes(subshapes);

    if (shapename == "triangle")
    {
        rawshapeptr = std::shared_ptr<rawtriangle>(new rawtriangle(physreg, subrawshapes, nummeshpts));
        return;
    }
    if (shapename == "quadrangle")
    {
        rawshapeptr = std::shared_ptr<rawquadrangle>(new rawquadrangle(physreg, subrawshapes, nummeshpts));
        return;
    }

    std::cout << "Error in 'shape' object: shape " << shapename << " does not accept this constructor or is unknown (try lower case)" << std::endl;
    abort();
}

shape::shape(std::string shapename, int physreg, std::vector<shape> subshapes)
{
    // Get the rawshape pointer from all shapes:
    std::vector<std::shared_ptr<rawshape>> subrawshapes = geotools::getrawshapes(subshapes);

    if (shapename == "triangle")
    {
        rawshapeptr = std::shared_ptr<rawtriangle>(new rawtriangle(physreg, subrawshapes));
        return;
    }
    if (shapename == "quadrangle")
    {
        rawshapeptr = std::shared_ptr<rawquadrangle>(new rawquadrangle(physreg, subrawshapes));
        return;
    }
    if (shapename == "union")
    {
        rawshapeptr = std::shared_ptr<rawunion>(new rawunion(physreg, subrawshapes));
        return;
    }

    std::cout << "Error in 'shape' object: shape " << shapename << " does not accept this constructor or is unknown (try lower case)" << std::endl;
    abort();
}

shape::shape(std::string shapename, int physreg, std::vector<double> centercoords, double radius, int nummeshpts)
{
    if (centercoords.size() != 3)
    {
        std::cout << "Error in 'shape' object: length of center point coordinate vector for shape " << shapename << " should be three" << std::endl;
        abort();
    }

    std::shared_ptr<rawshape> centerpoint = geotools::coordstopoints(centercoords)[0];

    if (shapename == "disk")
    {
        rawshapeptr = std::shared_ptr<rawdisk>(new rawdisk(physreg, centerpoint, radius, nummeshpts));
        return;
    }

    std::cout << "Error in 'shape' object: trying to call the disk constructor for a " << shapename << std::endl;
    abort();
}

shape::shape(std::string shapename, int physreg, shape centerpoint, double radius, int nummeshpts)
{
    if (shapename == "disk")
    {
        rawshapeptr = std::shared_ptr<rawdisk>(new rawdisk(physreg, centerpoint.getpointer(), radius, nummeshpts));
        return;
    }

    std::cout << "Error in 'shape' object: trying to call the disk constructor for a " << shapename << std::endl;
    abort();
}

void shape::setphysicalregion(int physreg)
{
    errornullpointer();
    rawshapeptr->setphysicalregion(physreg); 
}

void shape::move(expression u)
{
    errornullpointer();
    if (u.countcolumns() != 1 || u.countrows() > 3)
    {
        std::cout << "Error in 'shape' object: unexpected argument size in 'move' (expected a column vector up to length 3)" << std::endl;
        abort();
    }
    u = u.resize(3,1);
    rawshapeptr->move(mathop::compx(u), mathop::compy(u), mathop::compz(u));
}

void shape::shift(double shiftx, double shifty, double shiftz)
{
    errornullpointer();
    rawshapeptr->shift(shiftx, shifty, shiftz); 
}

void shape::scale(double scalex, double scaley, double scalez)
{
    errornullpointer();
    rawshapeptr->scale(scalex, scaley, scalez); 
}

void shape::rotate(double alphax, double alphay, double alphaz)
{
    errornullpointer();
    rawshapeptr->rotate(alphax, alphay, alphaz); 
}

shape shape::extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection)
{
    errornullpointer();

    if (numlayers < 2)
    {
        std::cout << "Error in 'shape' object: cannot extrude with " << numlayers << " node layers (at least two are needed)" << std::endl;
        abort();
    }
    if (extrudedirection.size() != 3)
    {
        std::cout << "Error in 'shape' object: extrude direction should be a vector of length three" << std::endl;
        abort();
    }
    
    myalgorithm::normvector(extrudedirection);
    
    return shape(rawshapeptr->duplicate()->extrude(physreg, height, numlayers, extrudedirection));
}

std::vector<shape> shape::extrude(std::vector<int> physreg, std::vector<double> height, std::vector<int> numlayers, std::vector<double> extrudedirection)
{
    if (physreg.size() != height.size() || physreg.size() != numlayers.size())
    {
        std::cout << "Error in 'shape' object: extrude vector arguments should have the same size" << std::endl;
        abort();
    }
 
    myalgorithm::normvector(extrudedirection);
    
    int num = physreg.size();
    std::vector<shape> output(num);

    std::vector<double> xyzshift = {0,0,0};
    for (int i = 0; i < num; i++)
    {
        shape curextr = shape(rawshapeptr->duplicate()).extrude(physreg[i], height[i], numlayers[i], extrudedirection);
        curextr.shift(xyzshift[0],xyzshift[1],xyzshift[2]);
        
        xyzshift[0] += height[i]*extrudedirection[0];
        xyzshift[1] += height[i]*extrudedirection[1];
        xyzshift[2] += height[i]*extrudedirection[2];
        
        output[i] = curextr;
    }
    
    return output;
}

shape shape::duplicate(void)
{
    errornullpointer();
    return shape(rawshapeptr->duplicate());
}

int shape::getdimension(void)
{
    errornullpointer();
    return rawshapeptr->getdimension(); 
}

std::vector<double> shape::getcoords(void)
{
    errornullpointer();
    return *(rawshapeptr->getcoords()); 
}

std::string shape::getname(void)
{
    errornullpointer();
    return rawshapeptr->getname(); 
}

std::vector<shape> shape::getsons(void)
{
    errornullpointer();
    return geotools::getshapes(rawshapeptr->getsons());
}

int shape::getphysicalregion(void)
{
    errornullpointer();
    return rawshapeptr->getphysicalregion(); 
}

std::shared_ptr<rawshape> shape::getpointer(void) 
{ 
    return rawshapeptr; 
}


