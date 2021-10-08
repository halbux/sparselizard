#include "rawshape.h"
#include "geotools.h"


void rawshape::move(expression xmove, expression ymove, expression zmove, bool recursively)
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

    std::vector<double> xmovevec = xmove.evaluate(xcoords,ycoords,zcoords);
    std::vector<double> ymovevec = ymove.evaluate(xcoords,ycoords,zcoords);
    std::vector<double> zmovevec = zmove.evaluate(xcoords,ycoords,zcoords);

    for (int i = 0; i < numnodes; i++)
    {
        mycoords->at(3*i+0) += xmovevec[i];
        mycoords->at(3*i+1) += ymovevec[i];
        mycoords->at(3*i+2) += zmovevec[i];
    }

    
    // Also move (only once) every sub shape:
    if (recursively)
    {
        std::vector<rawshape*> subshapes = geotools::unique(geotools::getpointers(getsubshapesrecursively()));
        for (int i = 0; i < subshapes.size(); i++)
            subshapes[i]->move(xmove, ymove, zmove, false);
    }
}

void rawshape::shift(double shiftx, double shifty, double shiftz, bool recursively)
{
    std::vector<double>* mycoords = getcoords();

    int numnodes = mycoords->size()/3;

    for (int i = 0; i < numnodes; i++)
    {
        mycoords->at(3*i+0) += shiftx;
        mycoords->at(3*i+1) += shifty;
        mycoords->at(3*i+2) += shiftz;
    }


    // Also shift (only once) every sub shape:
    if (recursively)
    {
        std::vector<rawshape*> subshapes = geotools::unique(geotools::getpointers(getsubshapesrecursively()));
        for (int i = 0; i < subshapes.size(); i++)
            subshapes[i]->shift(shiftx, shifty, shiftz, false);
    }
}

void rawshape::scale(double scalex, double scaley, double scalez, bool recursively)
{
    std::vector<double>* mycoords = getcoords();

    int numnodes = mycoords->size()/3;

    for (int i = 0; i < numnodes; i++)
    {
        mycoords->at(3*i+0) = scalex*mycoords->at(3*i+0);
        mycoords->at(3*i+1) = scaley*mycoords->at(3*i+1);
        mycoords->at(3*i+2) = scalez*mycoords->at(3*i+2);
    }


    // Also scale (only once) every sub shape:
    if (recursively)
    {
        std::vector<rawshape*> subshapes = geotools::unique(geotools::getpointers(getsubshapesrecursively()));
        for (int i = 0; i < subshapes.size(); i++)
            subshapes[i]->scale(scalex, scaley, scalez, false);
    }
}

void rawshape::rotate(double alphax, double alphay, double alphaz, bool recursively)
{
    // Get the coordinates of the shape to rotate:
    std::vector<double>* mycoords = getcoords();

    geotools::rotate(alphax, alphay, alphaz, mycoords);


    // Also rotate (only once) every sub shape:
    if (recursively)
    {
        std::vector<rawshape*> subshapes = geotools::unique(geotools::getpointers(getsubshapesrecursively()));
        for (int i = 0; i < subshapes.size(); i++)
            subshapes[i]->rotate(alphax, alphay, alphaz, false);
    }
}


std::shared_ptr<rawshape> rawshape::extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection)
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

std::vector<std::shared_ptr<rawshape>> rawshape::getsubshapesrecursively(void)
{
    std::vector<std::shared_ptr<rawshape>> output(countsubshapesrecursively());

    std::vector<std::shared_ptr<rawshape>> subshapes = getsubshapes();

    int index = 0;
    for (int i = 0; i < subshapes.size(); i++)
    {
        output[index] = subshapes[i];
        index++;

        std::vector<std::shared_ptr<rawshape>> toappend = subshapes[i]->getsubshapesrecursively();
        for (int j = 0; j < toappend.size(); j++)
            output[index+j] = toappend[j];
        index += toappend.size();
    }

    return output;
}

void rawshape::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    std::cout << "Error in 'rawshape' object: 'setsubshapes' has not been defined for this shape" << std::endl;
    abort(); 
}

void rawshape::setsubshapesrecursively(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    // Prepare to store the direct subshapes of this raw shape:
    int numberofdirectsubshapes = (getsubshapes()).size();
    std::vector<std::shared_ptr<rawshape>> directsubshapes(numberofdirectsubshapes);

    int subshapeindex = 0;
    for (int n = 0; n < numberofdirectsubshapes; n++)
    {
        directsubshapes[n] = subshapes[subshapeindex];
        subshapeindex++;
    
        // Count the number of subshapes recursively in the current direct subshape:
        int numberoflowerlevelsubshapes = directsubshapes[n]->countsubshapesrecursively();

        // Get all subshapes of lower levels:
        std::vector<std::shared_ptr<rawshape>> lowerlevelsubshapes(numberoflowerlevelsubshapes);
        for (int i = 0; i < numberoflowerlevelsubshapes; i++)
            lowerlevelsubshapes[i] = subshapes[subshapeindex+i];
        directsubshapes[n]->setsubshapesrecursively(lowerlevelsubshapes);
        
        subshapeindex += numberoflowerlevelsubshapes;
    }

    setsubshapes(directsubshapes);
}

int rawshape::countsubshapesrecursively(void)
{
    std::vector<std::shared_ptr<rawshape>> subshapes = getsubshapes();

    int numsubshapes = subshapes.size();

    for (int i = 0; i < subshapes.size(); i++)
        numsubshapes += subshapes[i]->countsubshapesrecursively();

    return numsubshapes;
}

void rawshape::replicatelinks(std::shared_ptr<rawshape> origin)
{
    // Get the rawshape pointers of all subshapes in the origin:
    std::vector<rawshape*> originsubshapes = geotools::getpointers(origin->getsubshapesrecursively());

    if (originsubshapes.size() <= 1)
        return;

    // Get all subshapes in this rawshape:
    std::vector<std::shared_ptr<rawshape>> subshapes = this->getsubshapesrecursively();

    // Copy the equality relations from the origin to 'subshapes':
    std::vector<int> reorderingvector;
    geotools::sortrawshapepointers(originsubshapes, reorderingvector);

    std::vector<std::shared_ptr<rawshape>> linkedsubshapes(reorderingvector.size());

    linkedsubshapes[reorderingvector[0]] = subshapes[reorderingvector[0]];
    for (int i = 1; i < reorderingvector.size(); i++)
    {
        if (originsubshapes[reorderingvector[i]] != originsubshapes[reorderingvector[i-1]])
            linkedsubshapes[reorderingvector[i]] = subshapes[reorderingvector[i]];
        else
            linkedsubshapes[reorderingvector[i]] = linkedsubshapes[reorderingvector[i-1]];
    }

    // Transfer linked subshapes to this rawshape:
    setsubshapesrecursively(linkedsubshapes);
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



