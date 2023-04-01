#include "geotools.h"


double geotools::acos(double arg)
{
    double pi = 3.1415926535897932384;

    // Extra margin because of the sensitivity of acos to noise at 1 and -1:
    double tol = 1e-10;
    if (arg >= 1.0-tol)
        return 0.0;
    if (arg <= -1.0+tol)
        return pi;

    return std::acos(arg);
}

std::vector<std::shared_ptr<rawshape>> geotools::coordstopoints(std::vector<double> coords)
{
    if (coords.size()%3 == 0)
    {
        int numpts = coords.size()/3;    

        std::vector<std::shared_ptr<rawshape>> pts(numpts);
        for (int i = 0; i < numpts; i++)
            pts[i] = std::shared_ptr<rawpoint>(new rawpoint(-1, {coords[3*i+0], coords[3*i+1], coords[3*i+2]}));

        return pts;
    }
    else
    {
        logs log;
        log.msg() << "Error in 'geotools' namespace: length of coordinate vector should be a multiple of three" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

double geotools::getdistance(std::vector<double> pt1coords, std::vector<double> pt2coords)
{
    return std::sqrt( std::pow(pt1coords[0]-pt2coords[0],2) + std::pow(pt1coords[1]-pt2coords[1],2) + std::pow(pt1coords[2]-pt2coords[2],2) );
}

double geotools::getdistance(int pt1, int pt2, std::vector<double>& coords)
{
    return std::sqrt( std::pow(coords[3*pt1+0]-coords[3*pt2+0],2) + std::pow(coords[3*pt1+1]-coords[3*pt2+1],2) + std::pow(coords[3*pt1+2]-coords[3*pt2+2],2) );
}

double geotools::getplanerotation(std::string xy, std::vector<double> p1, std::vector<double> p2, std::vector<double> p3)
{
    double pi = 3.1415926535897932384;

    // Roundoff noise threshold:
    double threshold = 1e-10;

    // Vector between two points:
    std::vector<double> v12 = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};
    std::vector<double> v13 = {p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]};
    double v12norm = std::sqrt(v12[0]*v12[0]+v12[1]*v12[1]+v12[2]*v12[2]);

    // Shortcut if all points are in a plane parallel to the xy plane.
    // This allows colinear points in 2D.
    if (std::abs(v12[2])/v12norm < threshold && std::abs(v13[2])/v12norm < threshold)
        return 0.0;

    // Normal [a b c] to the plane is the cross product:
    std::vector<double> vnormal = {v12[1]*v13[2]-v12[2]*v13[1], v12[2]*v13[0]-v12[0]*v13[2], v12[0]*v13[1]-v12[1]*v13[0]};
    double normalnorm = std::sqrt(vnormal[0]*vnormal[0]+vnormal[1]*vnormal[1]+vnormal[2]*vnormal[2]);

    // If the points are colinear:
    if (normalnorm/v12norm < threshold)
    {
        logs log;
        log.msg() << "Error in 'geotools' namespace: points provided are colinear (only allowed in 2D)";
        log.error();
    }

    // If the normal is perpendicular to the z axis (i.e. it has a 0 z component):
    if (std::abs(vnormal[2])/v12norm < threshold)
        return 90;

    // Plane equation is ax + by + cz = d.
    // Setting x or y to zero gives the angle we need:
    if (xy == "xrot")
        return -180.0/pi*std::atan(-vnormal[1]/vnormal[2]);
    if (xy == "yrot")
        return 180.0/pi*std::atan(-vnormal[0]/vnormal[2]);
        
    throw std::runtime_error(""); // fix return warning
}

void geotools::rotate(double alphax, double alphay, double alphaz, std::vector<double>* coords)
{
    int numberofnodes = coords->size()/3;

    // Convert input degrees to radians:
    double pi = 3.1415926535897932384;
    double ax = alphax*pi/180.0;
    double ay = alphay*pi/180.0;
    double az = alphaz*pi/180.0;
    
    // Define the rotation matrix R = Rz*Ry*Rx:
    double cx = std::cos(ax); double sx = std::sin(ax);
    double cy = std::cos(ay); double sy = std::sin(ay);
    double cz = std::cos(az); double sz = std::sin(az);
    
    double Rxx = cy*cz;
    double Rxy = cz*sx*sy - cx*sz;
    double Rxz = sx*sz + cx*cz*sy;
    double Ryx = cy*sz; 
    double Ryy = cx*cz + sx*sy*sz; 
    double Ryz = cx*sy*sz - cz*sx; 
    double Rzx = -sy; 
    double Rzy = cy*sx; 
    double Rzz = cx*cy; 

    // Compute R*[coordx; coordy; coordz]:
    for (int nodenumber = 0; nodenumber < numberofnodes; nodenumber++)
    {
        double xcoord = coords->at(3*nodenumber+0);
        double ycoord = coords->at(3*nodenumber+1);
        double zcoord = coords->at(3*nodenumber+2);
        
        coords->at(3*nodenumber+0) = Rxx * xcoord + Rxy * ycoord + Rxz * zcoord;
        coords->at(3*nodenumber+1) = Ryx * xcoord + Ryy * ycoord + Ryz * zcoord;
        coords->at(3*nodenumber+2) = Rzx * xcoord + Rzy * ycoord + Rzz * zcoord;
    }
}

std::vector<double> geotools::flipcoords(std::vector<double>& input)
{
    int numnodes = input.size()/3;

    std::vector<double> output(3*numnodes);

    for (int i = 0; i < numnodes; i++)
    {
        output[3*i+0] = input[3*(numnodes-1-i)+0];
        output[3*i+1] = input[3*(numnodes-1-i)+1];
        output[3*i+2] = input[3*(numnodes-1-i)+2];
    }
    
    return output;
}

std::vector<std::shared_ptr<rawshape>> geotools::orient(std::vector<std::shared_ptr<rawshape>> input)
{
    if (input.size() < 2)
        return input;

    // Get the coordinates of the points in the two first lines:
    std::vector<double> p1coord = *(input[0]->getsons()[0]->getcoords());
    std::vector<double> p2coord = *(input[0]->getsons()[1]->getcoords());
    std::vector<double> p3coord = *(input[1]->getsons()[0]->getcoords());
    std::vector<double> p4coord = *(input[1]->getsons()[1]->getcoords());

    // Find which pair between the two lines is the closest in distance:
    double d13 = getdistance(p1coord,p3coord);
    double d14 = getdistance(p1coord,p4coord);
    double d23 = getdistance(p2coord,p3coord);
    double d24 = getdistance(p2coord,p4coord);

    // This will tell which lines must be flipped:
    std::vector<bool> flipline(input.size(), false);

    if (d13 < d14 && d13 < d23 && d13 < d24)
        flipline[0] = true;
    if (d14 < d13 && d14 < d23 && d14 < d24)
        flipline[0] = flipline[1] = true;
    if (d24 < d13 && d24 < d23 && d24 < d14)
        flipline[1] = true;

    
    std::vector<double> prevnodecoord;
    if (flipline[1] == false)
        prevnodecoord = p4coord;
    else
        prevnodecoord = p3coord;

    for (int i = 2; i < input.size(); i++)
    {
        // Get the coordinates of the points in the current line:
        std::vector<double> p1c = *(input[i]->getsons()[0]->getcoords());
        std::vector<double> p2c = *(input[i]->getsons()[1]->getcoords());

        // Find which point is closest to the previous one:
        double d1 = getdistance(p1c,prevnodecoord);
        double d2 = getdistance(p2c,prevnodecoord);
            
        if (d1 < d2)
            prevnodecoord = p2c;
        else
        {
            flipline[i] = true;
            prevnodecoord = p1c;
        }
    }


    // Flip the lines that must be flipped:
    for (int i = 0; i < input.size(); i++)
    {
        if (flipline[i] == true)
            input[i]->flip();
    }

    return input;
}

std::vector< std::shared_ptr<rawshape> > geotools::getrawshapes(std::vector<shape> shapes)
{
    std::vector< std::shared_ptr<rawshape> > rawshapes(shapes.size());

    for (int i = 0; i < shapes.size(); i++)
    {
        std::shared_ptr<rawshape> currentrawshapeptr = shapes[i].getpointer();

        if (currentrawshapeptr != NULL)    
            rawshapes[i] = currentrawshapeptr;
        else
        {
            logs log;
            log.msg() << "Error in 'geotools' namespace: encountered an undefined shape (NULL rawshape pointer)" << std::endl;
            log.error();
        }
    }
    
    return rawshapes;
}

std::vector<shape> geotools::getshapes(std::vector< std::shared_ptr<rawshape> > rawshapes)
{
    std::vector<shape> shapes(rawshapes.size());

    for (int i = 0; i < rawshapes.size(); i++)
        shapes[i] = shape(rawshapes[i]);
    
    return shapes;
}

std::vector<rawshape*> geotools::getpointers(std::vector< std::shared_ptr<rawshape> > sharedptrs)
{
    std::vector<rawshape*> output(sharedptrs.size());
    for (int i = 0; i < output.size(); i++)
        output[i] = sharedptrs[i].get();

    return output;
}

std::vector< std::shared_ptr<rawshape> > geotools::flip(std::vector< std::shared_ptr<rawshape> > input)
{
    std::vector< std::shared_ptr<rawshape> > output(input.size());

    for (int i = 0; i < input.size(); i++)
        output[i] = input[input.size()-1-i];    

    return output;
}

std::vector<rawshape*> geotools::unique(std::vector<rawshape*> ptrs)
{
    std::sort(ptrs.begin(), ptrs.end());
    ptrs.erase( std::unique(ptrs.begin(), ptrs.end()), ptrs.end());

    return ptrs;
}

std::vector<std::shared_ptr<rawshape>> geotools::duplicate(std::vector<std::shared_ptr<rawshape>> input)
{
    std::vector<std::shared_ptr<rawshape>> output(input.size());

    for (int i = 0; i < input.size(); i++)
        output[i] = input[i]->duplicate();

    return output;
}

std::vector<std::shared_ptr<rawshape>> geotools::concatenate(std::vector<std::vector<std::shared_ptr<rawshape>>> input)
{
    // Compute the total number of rawshapes:
    int numrawshapes = 0;
    for (int i = 0; i < input.size(); i++)
        numrawshapes += input[i].size();

    std::vector<std::shared_ptr<rawshape>> output(numrawshapes);

    int index = 0;
    for (int i = 0; i < input.size(); i++)
    {
        for (int j = 0; j < input[i].size(); j++)
        {
            output[index] = input[i][j];
            index++;
        }
    }

    return output;
}

void geotools::sortrawshapepointers(std::vector<rawshape*>& tosort, std::vector<int>& reorderingvector)
{
    if (reorderingvector.size() != tosort.size())
        reorderingvector.resize(tosort.size());
    
    // Set 'reorderingvector' to [0 1 2 ...]:
    std::iota(reorderingvector.begin(), reorderingvector.end(), 0);
    // Sort 'reorderingvector' according to 'tosort':
    // The < operator is overloaded by a lambda function.
    std::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
        { 
            if (tosort[elem1] < tosort[elem2])
                return true;
            if (tosort[elem1] > tosort[elem2])
                return false;
            // For identical entries make a COHERENT decision for a stable sorting.
            return elem1 < elem2;
        });
}

std::vector<double> geotools::appendcoords(std::vector<std::shared_ptr<rawshape>> rawshapes)
{
    std::vector< std::vector<double>* > coordsptrs(rawshapes.size());

    for (int i = 0; i < rawshapes.size(); i++)
        coordsptrs[i] = rawshapes[i]->getcoords();
    

    // First calculate the total size:
    int totallen = 0;
    for (int i = 0; i < coordsptrs.size(); i++)
        totallen += coordsptrs[i]->size();

    std::vector<double> output(totallen);

    int index = 0;
    for (int i = 0; i < coordsptrs.size(); i++)
    {
        for (int j = 0; j < coordsptrs[i]->size(); j++)
        {
            output[index] = coordsptrs[i]->at(j);
            index++;
        }
    }
    return output;
}

std::vector<std::vector<int>> geotools::appendelems(std::vector<std::shared_ptr<rawshape>> rawshapes)
{
    std::vector< std::vector<std::vector<int>>* > elemsptrs(rawshapes.size());

    for (int i = 0; i < rawshapes.size(); i++)
        elemsptrs[i] = rawshapes[i]->getelems();
    

    std::vector<std::vector<int>> output(8);

    // Loop on all elements:
    for (int e = 0; e < 8; e++)
    {
        // First calculate the total size:
        int totallen = 0;
        for (int i = 0; i < elemsptrs.size(); i++)
            totallen += elemsptrs[i]->at(e).size();

        output[e] = std::vector<int>(totallen);

        int index = 0, nodenumshift = 0;
        for (int i = 0; i < elemsptrs.size(); i++)
        {
            for (int j = 0; j < elemsptrs[i]->at(e).size(); j++)
            {
                output[e][index] = elemsptrs[i]->at(e)[j] + nodenumshift;
                index++;
            }
            nodenumshift += rawshapes[i]->getcoords()->size()/3;
        }
    }
    return output;
}

int geotools::getcurvatureorder(std::vector< std::shared_ptr<rawshape> > rawshapes)
{
    if (rawshapes.size() == 0)
        return -1;
        
    int co = rawshapes[0]->getcurvatureorder();
    
    for (int i = 1; i < rawshapes.size(); i++)
    {
        int curco = rawshapes[i]->getcurvatureorder();
        if (curco != co)
        {
            logs log;
            log.msg() << "Error in 'geotools' namespace: expected a unique curvature order for all shapes provided" << std::endl;
            log.error();
        }
    }
    
    return co;
}

