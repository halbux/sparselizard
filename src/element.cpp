#include "element.h"
#include "geotools.h"
#include "lagrangeformfunction.h"


namespace sl {

element::element(std::string elementname)
{
    curvedtypenumber = -1;
    // 'switch' does not work on std::string:
    if (elementname == "point")
        curvedtypenumber = 0;
    if (elementname == "line")
        curvedtypenumber = 1;
    if (elementname == "triangle")
        curvedtypenumber = 2;
    if (elementname == "quadrangle")
        curvedtypenumber = 3;
    if (elementname == "tetrahedron")
        curvedtypenumber = 4;
    if (elementname == "hexahedron")
        curvedtypenumber = 5;
    if (elementname == "prism")
        curvedtypenumber = 6;
    if (elementname == "pyramid")
        curvedtypenumber = 7;
    if (curvedtypenumber == -1)
    {
        logs log;
        log.msg() << "Error in 'element' object: trying to use undefined element name: " << elementname << std::endl << "Make sure everything is lower case" << std::endl;
        log.error();
    }
}

element::element(int number)
{
    if (number < 0)
    {
        logs log;
        log.msg() << "Error in 'element' object: cannot define negative element type number " << number << std::endl;
        log.error();
    }
    curvedtypenumber = number;
}

element::element(int number, int curvatureorder)
{
    if (number < 0)
    {
        logs log;
        log.msg() << "Error in 'element' object: can not define a negative element type number" << std::endl;
        log.msg() << "Element type number is " << number << " with " << curvatureorder << " curvature order" << std::endl;
        log.error();
    }
    if (curvatureorder <= 0)
    {
        logs log;
        log.msg() << "Error in 'element' object: can not define a negative or 0 curvature order" << std::endl;
        log.msg() << "Element type number is " << number << " with " << curvatureorder << " curvature order" << std::endl;
        log.error();
    }
    // The point element can only have number 0:
    if (number == 0)
        curvedtypenumber = 0;
    else
        curvedtypenumber = 7*(curvatureorder-1)+number;
}

void element::setnodes(std::vector<int>& nodelist)
{
    if (curvedtypenumber == -1)
    {
        logs log;
        log.msg() << "Error: element type has not been defined yet" << std::endl;
        log.error();
    }
    if (nodelist.size() != countcurvednodes())
    {
        logs log;
        log.msg() << "Error: trying to define an order " << getcurvatureorder() << " " << gettypename() << " with " << nodelist.size() << " nodes. There should be " << countcurvednodes() << std::endl;
        log.error();
    }
    curvednodelist = nodelist;
}

std::vector<int> element::getnodes(void)
{
    return curvednodelist;
}

std::string element::gettypename(void)
{
    switch (gettypenumber())
    {
        case 0:
            return "point";
        case 1:
            return "line";
        case 2:
            return "triangle";
        case 3:
            return "quadrangle";
        case 4:
            return "tetrahedron";
        case 5:
            return "hexahedron";
        case 6:
            return "prism";
        case 7:
            return "pyramid";
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::string element::gettypenameconjugation(int numberofelements)
{
    if (numberofelements < 2)
        return gettypename();
    else
    {
        switch (gettypenumber())
        {
            case 0:
                return "points";
            case 1:
                return "lines";
            case 2:
                return "triangles";
            case 3:
                return "quadrangles";
            case 4:
                return "tetrahedra";
            case 5:
                return "hexahedra";
            case 6:
                return "prisms";
            case 7:
                return "pyramids";
        }
    }
    
    throw std::runtime_error(""); // fix return warning
}

bool element::iscurved(void)
{
    return (curvedtypenumber > 7);
}

int element::getcurvatureorder(void)
{
    if (curvedtypenumber%7 == 0 && curvedtypenumber != 0)
        return curvedtypenumber/7;
    else
        return (curvedtypenumber - curvedtypenumber%7)/7 + 1;        
}

int element::gettypenumber(void)
{
    return (curvedtypenumber - 7*(getcurvatureorder() - 1));
}

int element::getcurvedtypenumber(void)
{
    return curvedtypenumber;
}

int element::countcurvednodes(void)
{
    int order = getcurvatureorder();
    int curvednumberofnodes;

    switch (gettypenumber())
    {
        // Point:
        case 0:
            curvednumberofnodes = 1;
            break;
        // Line:
        case 1:
            curvednumberofnodes = order + 1;
            break;
        // Triangle:
        case 2:
            // (#on quad + #on diagonal) / 2:
            curvednumberofnodes = ( std::pow(order+1,2) + (order+1) )/2;
            break;
        // Quadrangle:
        case 3:
            curvednumberofnodes = std::pow(order + 1,2);
            break;
        // Tetrahedron:
        case 4:
            curvednumberofnodes = ( (order+1)*(order+2)*(order+3) )/6;
            break;
        // Hexahedron:
        case 5:
            curvednumberofnodes = std::pow(order + 1,3);
            break;
        // Prism:
        case 6:
            curvednumberofnodes = ( (order + 1)*( std::pow(order+1,2) + (order+1) ) )/2;
            break;
        // Pyramid:
        case 7:
            // sum of all parallel quadrangles:
            curvednumberofnodes = 0;
            for (int i = 0; i < order + 1; i++)
                curvednumberofnodes = curvednumberofnodes + std::pow(i + 1,2);
            break;
    }

    return curvednumberofnodes;
}

int element::getelementdimension(void)
{
    int straighttypenumber = gettypenumber();
    
    if (straighttypenumber == 0)
        return 0;
    if (straighttypenumber == 1)
        return 1;
    if (straighttypenumber == 2 || straighttypenumber == 3)
        return 2;
    if (straighttypenumber > 3)
        return 3;
        
    throw std::runtime_error(""); // fix return warning
}


int element::counttype(int typenum)
{
    switch (typenum)
    {
        // Point:
        case 0:
        {
            switch (gettypenumber())
            {
                // Point:
                case 0:
                    return 1;
                // Line:
                case 1:
                    return 2;
                // Triangle:
                case 2:
                    return 3;
                // Quadrangle:
                case 3:
                    return 4;
                // Tetrahedron:
                case 4:
                    return 4;
                // Hexahedron:
                case 5:
                    return 8;
                // Prism:
                case 6:
                    return 6;
                // Pyramid:
                case 7:
                    return 5;
            }
        }
        // Line:
        case 1:
        {
            switch (gettypenumber())
            {
                // Point:
                case 0:
                    return 0;
                // Line:
                case 1:
                    return 1;
                // Triangle:
                case 2:
                    return 3;
                // Quadrangle:
                case 3:
                    return 4;
                // Tetrahedron:
                case 4:
                    return 6;
                // Hexahedron:
                case 5:
                    return 12;
                // Prism:
                case 6:
                    return 9;
                // Pyramid:
                case 7:
                    return 8;
            }
        }
        // Triangle:
        case 2:
        {
            switch (gettypenumber())
            {
                // Point:
                case 0:
                    return 0;
                // Line:
                case 1:
                    return 0;
                // Triangle:
                case 2:
                    return 1;
                // Quadrangle:
                case 3:
                    return 0;
                // Tetrahedron:
                case 4:
                    return 4;
                // Hexahedron:
                case 5:
                    return 0;
                // Prism:
                case 6:
                    return 2;
                // Pyramid:
                case 7:
                    return 4;
            }
        }
        // Quadrangle:
        case 3:
        {
            switch (gettypenumber())
            {
                // Point:
                case 0:
                    return 0;
                // Line:
                case 1:
                    return 0;
                // Triangle:
                case 2:
                    return 0;
                // Quadrangle:
                case 3:
                    return 1;
                // Tetrahedron:
                case 4:
                    return 0;
                // Hexahedron:
                case 5:
                    return 6;
                // Prism:
                case 6:
                    return 3;
                // Pyramid:
                case 7:
                    return 1;
            }
        }
        // Tetrahedron:
        case 4:
        {
            if (gettypenumber() == 4)
                return 1;
            else
                return 0;
        }
        // Hexahedron:
        case 5:
        {
            if (gettypenumber() == 5)
                return 1;
            else
                return 0;
        }
        // Prism:
        case 6:
        {
            if (gettypenumber() == 6)
                return 1;
            else
                return 0;
        }
        // Pyramid:
        case 7:
        {
            if (gettypenumber() == 7)
                return 1;
            else
                return 0;
        }
    }
    
    throw std::runtime_error(""); // fix return warning
}

int element::countdim(int dim)
{
    switch (dim)
    {
        case 0:
            return counttype(0);
        case 1:
            return counttype(1);
        case 2:
            return counttype(2)+counttype(3);
        case 3:
            return counttype(4)+counttype(5)+counttype(6)+counttype(7);
    }
    
    throw std::runtime_error(""); // fix return warning
}

int element::countnodes(void)
{
    return counttype(0);
}

int element::countedges(void)
{
    return counttype(1);
}

int element::countfaces(void)
{
    return countdim(2);
}

int element::counttriangularfaces(void)
{
    return counttype(2);
}

int element::countquadrangularfaces(void)
{
    return counttype(3);
}

int element::countvolumes(void)
{
    return countdim(3);
}


bool element::isinsideelement(double ki, double eta, double phi)
{
    double roundoffnoise = 1e-10;

    switch (gettypenumber())
    {
        // Point:
        case 0:
            return (std::abs(ki) < roundoffnoise && std::abs(eta) < roundoffnoise && std::abs(phi) < roundoffnoise);
        // Line:
        case 1:
            return (std::abs(ki) < 1+roundoffnoise && std::abs(eta) < roundoffnoise && std::abs(phi) < roundoffnoise);
        // Triangle:
        case 2:
            return (ki+eta < 1+roundoffnoise && ki > -roundoffnoise && eta > -roundoffnoise && std::abs(phi) < roundoffnoise);
        // Quadrangle:
        case 3:
            return (std::abs(ki) < 1+roundoffnoise && std::abs(eta) < 1+roundoffnoise && std::abs(phi) < roundoffnoise);
        // Tetrahedron:
        case 4:
            return (ki+eta+phi < 1+roundoffnoise && ki > -roundoffnoise && eta > -roundoffnoise && phi > -roundoffnoise);
        // Hexahedron:
        case 5:
            return (std::abs(ki) < 1+roundoffnoise && std::abs(eta) < 1+roundoffnoise && std::abs(phi) < 1+roundoffnoise);
        // Prism:
        case 6:
            return (ki+eta < 1+roundoffnoise && ki > -roundoffnoise && eta > -roundoffnoise && std::abs(phi) < 1+roundoffnoise);
        // Pyramid:
        case 7:
            return (std::abs(ki) < 1-phi+roundoffnoise && std::abs(eta) < 1-phi+roundoffnoise && phi > -roundoffnoise && phi < 1+roundoffnoise);
    }
    
    throw std::runtime_error(""); // fix return warning
}

void element::isinsideelement(std::vector<double>& coords, std::vector<double>& cc, std::vector<bool>& isinside, double roundoffnoise)
{
    int numcoords = coords.size()/3;
    isinside = std::vector<bool>(numcoords, true);

    int dim = getelementdimension();
    if (dim == 1)
    {   
        double deltax = std::abs(cc[3]-cc[0]) + roundoffnoise;
        for (int i = 0; i < numcoords; i++)
            isinside[i] = (std::abs(coords[3*i+0]-cc[0]) < deltax && std::abs(coords[3*i+0]-cc[3]) < deltax);
    }
    if (dim == 2)
    {   
        // Normal to each edge and pointing inside element:
        int ne = countedges();
        std::vector<double> normals(2*ne);
        std::vector<double> rofn(ne);
        // Calculate element orientation:
        std::vector<double> va = {cc[3*1+0]-cc[0], cc[3*1+1]-cc[1]};
        std::vector<double> vb = {cc[3*(ne-1)+0]-cc[0], cc[3*(ne-1)+1]-cc[1]};
        bool isanticw = ((va[0]*vb[1]-va[1]*vb[0]) > 0); // z comp of cross product
        
        for (int i = 0; i < ne; i++)
        {
            int ni = (i+1)%ne;
            double dx = cc[3*ni+0]-cc[3*i+0];
            double dy = cc[3*ni+1]-cc[3*i+1];
            if (not(isanticw))
            {
                dx = -dx;
                dy = -dy;
            }
            normals[2*i+0] = -dy;
            normals[2*i+1] = dx;
            rofn[i] = roundoffnoise*(std::abs(dx)+std::abs(dy));
        }
        for (int i = 0; i < numcoords; i++)
        {
            for (int j = 0; j < ne; j++)
            {
                double dx = coords[3*i+0]-cc[3*j+0];
                double dy = coords[3*i+1]-cc[3*j+1];
                // Scalar product must be positive:
                double sp = normals[2*j+0]*dx + normals[2*j+1]*dy;
                if (sp < -rofn[j])
                {
                    isinside[i] = false;
                    break;
                }
            }
        }
    }
    if (dim == 3)
    {   
        // Normal to each face and pointing outside element:
        int ntf = counttriangularfaces();
        int nf = countfaces();
        std::vector<int> facedef = getfacesdefinitionsbasedonnodes();
        
        std::vector<double> normals(3*nf);
        std::vector<double> rofn(nf);
        std::vector<double> firstnodeinface(3*nf);
        // Triangular faces come first:
        int fdi = 0;
        for (int i = 0; i < nf; i++)
        {
            int no = facedef[fdi+0];
            int na = facedef[fdi+1];
            int nb;
            if (i < ntf)
                nb = facedef[fdi+2];
            else
                nb = facedef[fdi+3];
            
            // Cross-product:
            std::vector<double> a = {cc[3*na+0]-cc[3*no+0], cc[3*na+1]-cc[3*no+1], cc[3*na+2]-cc[3*no+2]};
            std::vector<double> b = {cc[3*nb+0]-cc[3*no+0], cc[3*nb+1]-cc[3*no+1], cc[3*nb+2]-cc[3*no+2]};
            normals[3*i+0] = a[1]*b[2]-a[2]*b[1];
            normals[3*i+1] = a[2]*b[0]-a[0]*b[2];
            normals[3*i+2] = a[0]*b[1]-a[1]*b[0];
            rofn[i] = roundoffnoise * (std::abs(a[0])+std::abs(a[1])+std::abs(a[2])) * (std::abs(b[0])+std::abs(b[1])+std::abs(b[2]));
            
            firstnodeinface[3*i+0] = cc[3*no+0];
            firstnodeinface[3*i+1] = cc[3*no+1];
            firstnodeinface[3*i+2] = cc[3*no+2];
            
            if (i < ntf)
                fdi += 3;
            else
                fdi += 4;
        }
        for (int i = 0; i < numcoords; i++)
        {
            for (int j = 0; j < nf; j++)
            {
                double dx = coords[3*i+0]-firstnodeinface[3*j+0];
                double dy = coords[3*i+1]-firstnodeinface[3*j+1];
                double dz = coords[3*i+2]-firstnodeinface[3*j+2];
                // Scalar product must be negative:
                double sp = normals[3*j+0]*dx + normals[3*j+1]*dy + normals[3*j+2]*dz;
                if (sp > rofn[j])
                {
                    isinside[i] = false;
                    break;
                }
            }
        }
    }
}

void element::atnode(std::vector<double>& refcoords, std::vector<int>& nodenums, double roundoffnoise)
{
    int numrefs = refcoords.size()/3;
    nodenums = std::vector<int>(numrefs, -1);

    switch (gettypenumber())
    {
        case 1:
        {
            for (int i = 0; i < numrefs; i++)
            {
                if (std::abs(refcoords[3*i+0] + 1.0) < roundoffnoise)
                    nodenums[i] = 0;
                if (std::abs(refcoords[3*i+0] - 1.0) < roundoffnoise)
                    nodenums[i] = 1; 
            }
            break;
        }
        case 2:
        {
            for (int i = 0; i < numrefs; i++)
            {
                if (std::abs(refcoords[3*i+0] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+1] + 0.0) < roundoffnoise)
                    nodenums[i] = 0;
                if (std::abs(refcoords[3*i+0] - 1.0) < roundoffnoise && std::abs(refcoords[3*i+1] + 0.0) < roundoffnoise)
                    nodenums[i] = 1; 
                if (std::abs(refcoords[3*i+0] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+1] - 1.0) < roundoffnoise)
                    nodenums[i] = 2;
            }
            break;
        }
        case 3:
        {
            for (int i = 0; i < numrefs; i++)
            {
                if (std::abs(refcoords[3*i+0] + 1.0) < roundoffnoise && std::abs(refcoords[3*i+1] + 1.0) < roundoffnoise)
                    nodenums[i] = 0;
                if (std::abs(refcoords[3*i+0] - 1.0) < roundoffnoise && std::abs(refcoords[3*i+1] + 1.0) < roundoffnoise)
                    nodenums[i] = 1; 
                if (std::abs(refcoords[3*i+0] - 1.0) < roundoffnoise && std::abs(refcoords[3*i+1] - 1.0) < roundoffnoise)
                    nodenums[i] = 2;
                if (std::abs(refcoords[3*i+0] + 1.0) < roundoffnoise && std::abs(refcoords[3*i+1] - 1.0) < roundoffnoise)
                    nodenums[i] = 3;
            }
            break;
        }
        case 4:
        {
            for (int i = 0; i < numrefs; i++)
            {
                if (std::abs(refcoords[3*i+0] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+1] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+2] + 0.0) < roundoffnoise)
                    nodenums[i] = 0;
                if (std::abs(refcoords[3*i+0] - 1.0) < roundoffnoise && std::abs(refcoords[3*i+1] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+2] + 0.0) < roundoffnoise)
                    nodenums[i] = 1; 
                if (std::abs(refcoords[3*i+0] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+1] - 1.0) < roundoffnoise && std::abs(refcoords[3*i+2] + 0.0) < roundoffnoise)
                    nodenums[i] = 2;
                if (std::abs(refcoords[3*i+0] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+1] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+2] - 1.0) < roundoffnoise)
                    nodenums[i] = 3;
            }
            break;
        }
        // NOT DEFINED FOR HEXAHEDRA, PRISMS AND PYRAMIDS YET
    }
}

void element::atedge(std::vector<double>& refcoords, std::vector<int>& edgenums, double roundoffnoise)
{
    int numrefs = refcoords.size()/3;
    edgenums = std::vector<int>(numrefs, -1);

    switch (gettypenumber())
    {
        case 2:
        {
            for (int i = 0; i < numrefs; i++)
            {
                if (std::abs(refcoords[3*i+1] + 0.0) < roundoffnoise)
                    edgenums[i] = 0;
                if (std::abs(refcoords[3*i+0] + refcoords[3*i+1] - 1.0) < roundoffnoise)
                    edgenums[i] = 1; 
                if (std::abs(refcoords[3*i+0] + 0.0) < roundoffnoise)
                    edgenums[i] = 2;
            }
            break;
        }
        case 3:
        {
            for (int i = 0; i < numrefs; i++)
            {
                if (std::abs(refcoords[3*i+1] + 1.0) < roundoffnoise)
                    edgenums[i] = 0;
                if (std::abs(refcoords[3*i+0] - 1.0) < roundoffnoise)
                    edgenums[i] = 1; 
                if (std::abs(refcoords[3*i+1] - 1.0) < roundoffnoise)
                    edgenums[i] = 2;
                if (std::abs(refcoords[3*i+0] + 1.0) < roundoffnoise)
                    edgenums[i] = 3;
            }
            break;
        }
        case 4:
        {
            for (int i = 0; i < numrefs; i++)
            {
                if (std::abs(refcoords[3*i+1] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+2] + 0.0) < roundoffnoise)
                    edgenums[i] = 0;
                if (std::abs(refcoords[3*i+0] + refcoords[3*i+1] - 1.0) < roundoffnoise)
                    edgenums[i] = 1; 
                if (std::abs(refcoords[3*i+0] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+2] + 0.0) < roundoffnoise)
                    edgenums[i] = 2;
                if (std::abs(refcoords[3*i+0] + 0.0) < roundoffnoise && std::abs(refcoords[3*i+1] + 0.0) < roundoffnoise)
                    edgenums[i] = 3;
                if (std::abs(refcoords[3*i+1] + refcoords[3*i+2] - 1.0) < roundoffnoise)
                    edgenums[i] = 4;
                if (std::abs(refcoords[3*i+0] + refcoords[3*i+2] - 1.0) < roundoffnoise)
                    edgenums[i] = 5;
            }
            break;
        }
        // NOT DEFINED FOR HEXAHEDRA, PRISMS AND PYRAMIDS YET
    }
}

void element::atface(std::vector<double>& refcoords, std::vector<int>& facenums, double roundoffnoise)
{
    int numrefs = refcoords.size()/3;
    facenums = std::vector<int>(numrefs, -1);

    switch (gettypenumber())
    {
        case 4:
        {
            for (int i = 0; i < numrefs; i++)
            {
                if (std::abs(refcoords[3*i+2] + 0.0) < roundoffnoise)
                    facenums[i] = 0;
                if (std::abs(refcoords[3*i+1] + 0.0) < roundoffnoise)
                    facenums[i] = 1; 
                if (std::abs(refcoords[3*i+0] + 0.0) < roundoffnoise)
                    facenums[i] = 2;
                if (std::abs(refcoords[3*i+0] + refcoords[3*i+1] + refcoords[3*i+2] - 1.0) < roundoffnoise)
                    facenums[i] = 3;
            }
            break;
        }
        // NOT DEFINED FOR HEXAHEDRA, PRISMS AND PYRAMIDS YET
    }
}

std::vector<double> element::getedgebarycenter(std::vector<double>& nc)
{
    switch (curvedtypenumber)
    {
        case 1:
            return {0.5*(nc[3*0+0]+nc[3*1+0]), 0.5*(nc[3*0+1]+nc[3*1+1]), 0.5*(nc[3*0+2]+nc[3*1+2])};
        case 2:
            return {0.5*(nc[3*0+0]+nc[3*1+0]), 0.5*(nc[3*0+1]+nc[3*1+1]), 0.5*(nc[3*0+2]+nc[3*1+2]),
                    0.5*(nc[3*1+0]+nc[3*2+0]), 0.5*(nc[3*1+1]+nc[3*2+1]), 0.5*(nc[3*1+2]+nc[3*2+2]),
                    0.5*(nc[3*2+0]+nc[3*0+0]), 0.5*(nc[3*2+1]+nc[3*0+1]), 0.5*(nc[3*2+2]+nc[3*0+2])};
        case 3:
            return {0.5*(nc[3*0+0]+nc[3*1+0]), 0.5*(nc[3*0+1]+nc[3*1+1]), 0.5*(nc[3*0+2]+nc[3*1+2]),
                    0.5*(nc[3*1+0]+nc[3*2+0]), 0.5*(nc[3*1+1]+nc[3*2+1]), 0.5*(nc[3*1+2]+nc[3*2+2]),
                    0.5*(nc[3*2+0]+nc[3*3+0]), 0.5*(nc[3*2+1]+nc[3*3+1]), 0.5*(nc[3*2+2]+nc[3*3+2]),
                    0.5*(nc[3*3+0]+nc[3*0+0]), 0.5*(nc[3*3+1]+nc[3*0+1]), 0.5*(nc[3*3+2]+nc[3*0+2])};
        case 4:
            return {0.5*(nc[3*0+0]+nc[3*1+0]), 0.5*(nc[3*0+1]+nc[3*1+1]), 0.5*(nc[3*0+2]+nc[3*1+2]),
                    0.5*(nc[3*1+0]+nc[3*2+0]), 0.5*(nc[3*1+1]+nc[3*2+1]), 0.5*(nc[3*1+2]+nc[3*2+2]),
                    0.5*(nc[3*2+0]+nc[3*0+0]), 0.5*(nc[3*2+1]+nc[3*0+1]), 0.5*(nc[3*2+2]+nc[3*0+2]),
                    0.5*(nc[3*3+0]+nc[3*0+0]), 0.5*(nc[3*3+1]+nc[3*0+1]), 0.5*(nc[3*3+2]+nc[3*0+2]),
                    0.5*(nc[3*3+0]+nc[3*2+0]), 0.5*(nc[3*3+1]+nc[3*2+1]), 0.5*(nc[3*3+2]+nc[3*2+2]),
                    0.5*(nc[3*3+0]+nc[3*1+0]), 0.5*(nc[3*3+1]+nc[3*1+1]), 0.5*(nc[3*3+2]+nc[3*1+2])};
        // NOT DEFINED FOR HEXAHEDRA, PRISMS AND PYRAMIDS YET
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<double> element::getfacebarycenter(std::vector<double>& nc)
{
    switch (curvedtypenumber)
    {
        case 2:
            return {(nc[3*0+0]+nc[3*1+0]+nc[3*2+0])/3.0, (nc[3*0+1]+nc[3*1+1]+nc[3*2+1])/3.0, (nc[3*0+2]+nc[3*1+2]+nc[3*2+2])/3.0};
        case 3:
            return {(nc[3*0+0]+nc[3*1+0]+nc[3*2+0]+nc[3*3+0])/4.0, (nc[3*0+1]+nc[3*1+1]+nc[3*2+1]+nc[3*3+1])/4.0, (nc[3*0+2]+nc[3*1+2]+nc[3*2+2]+nc[3*3+2])/4.0};
        case 4:
            return {(nc[3*0+0]+nc[3*2+0]+nc[3*1+0])/3.0, (nc[3*0+1]+nc[3*2+1]+nc[3*1+1])/3.0, (nc[3*0+2]+nc[3*2+2]+nc[3*1+2])/3.0,
                    (nc[3*0+0]+nc[3*1+0]+nc[3*3+0])/3.0, (nc[3*0+1]+nc[3*1+1]+nc[3*3+1])/3.0, (nc[3*0+2]+nc[3*1+2]+nc[3*3+2])/3.0,
                    (nc[3*0+0]+nc[3*3+0]+nc[3*2+0])/3.0, (nc[3*0+1]+nc[3*3+1]+nc[3*2+1])/3.0, (nc[3*0+2]+nc[3*3+2]+nc[3*2+2])/3.0,
                    (nc[3*3+0]+nc[3*1+0]+nc[3*2+0])/3.0, (nc[3*3+1]+nc[3*1+1]+nc[3*2+1])/3.0, (nc[3*3+2]+nc[3*1+2]+nc[3*2+2])/3.0};
        // NOT DEFINED FOR HEXAHEDRA, PRISMS AND PYRAMIDS YET
    }
    
    throw std::runtime_error(""); // fix return warning
}

double element::measurereferenceelement(void)
{
    switch (gettypenumber())
    {
        // Point:
        case 0:
            return 1.0;
        // Line:
        case 1:
            return 2.0;
        // Triangle:
        case 2:
            return 0.5;
        // Quadrangle:
        case 3:
            return 4.0;
        // Tetrahedron:
        case 4:
            return 1.0/6.0;
        // Hexahedron:
        case 5:
            return 8.0;
        // Prism:
        case 6:
            return 1.0;
        // Pyramid:
        case 7:
            return 4.0/3.0;
    }
    
    throw std::runtime_error(""); // fix return warning
}

bool element::istriangularface(int facenum)
{
    switch (gettypenumber())
    {
        // Point:
        case 0:
            return false;
        // Line:
        case 1:
            return false;
        // Triangle:
        case 2:
            return true;
        // Quadrangle:
        case 3:
            return false;
        // Tetrahedron:
        case 4:
            return true;
        // Hexahedron:
        case 5:
            return false;
        // Prism:
        case 6:
            if (facenum < 2)
                return true;
            else
                return false;
        // Pyramid:
        case 7:
            if (facenum < 4)
                return true;
            else
                return false;
    }
    
    throw std::runtime_error(""); // fix return warning
}

bool element::ishorizontaledge(int edgenum)
{
    if (edgenum == 2 || edgenum == 4 || edgenum == 5)
        return false;
    else
        return true;
}

std::vector<int> element::getnodesinline(int lineindex)
{
    int order = getcurvatureorder();
    int numberofnodesinline = order + 1;
    std::vector<int> nodesinline(numberofnodesinline);
    int numberofinteriorlinenodes = numberofnodesinline - 2;
    
    std::vector<int> cornernodesinalledges = getedgesdefinitionsbasedonnodes();
    
    // The first nodes in the node list of the line are the corner nodes, followed by the interior nodes.
    // All interior nodes are listed consecutively in the node list of the element object from which we extract the lines.
    // Moreover, since we follow the edges in the correct direction the nodes appear in the correct order and the
    // order must thus not be reversed.
    nodesinline[0] = curvednodelist[cornernodesinalledges[2*lineindex+0]];
    nodesinline[1] = curvednodelist[cornernodesinalledges[2*lineindex+1]];
    
    int index = 2;
    for (int i = countnodes() + lineindex*numberofinteriorlinenodes; i < countnodes()+(lineindex+1)*numberofinteriorlinenodes; i++)
    {
        nodesinline[index] = curvednodelist[i];
        index = index + 1;
    }

    return nodesinline;
}

std::vector<int> element::getnodesintriangle(int triangleindex)
{
    return getnodesinsurface(triangleindex, true, false);
}

std::vector<int> element::getnodesinquadrangle(int quadrangleindex)
{
    return getnodesinsurface(quadrangleindex, false, true);
}

// 'getnodesinsurface' is only to be used by 'getnodesintriangle' and 'getnodesinquadrangle'.
// In 'getnodesinsurface' the nodes of the element object that correspond to its
// 'surfaceindex'th face are extracted and returned in an appropriate order
// to define the extracted surface.
// The general order for the nodes in curved lines, surfaces and volumes is:
//
// 1. The corner nodes in the right order
// 2. The interior nodes of edge 1 followed in the edge orientation then edge 2, edge 3,...
// 3. The interior nodes of surface 1 then surface 2,... (triangles before quadrangles)
// 4. The interior nodes of the volume
std::vector<int> element::getnodesinsurface(int surfaceindex, bool faceistriangle, bool faceisquadrangle)
{
    int order = getcurvatureorder();
    int numberofnodesinline = order + 1;
    int numberofinteriorlinenodes = numberofnodesinline - 2;
    int numberofnodesinsurface;
    int numberofcornernodesinsurface;
    int numberofedgesinsurface;
    int offset;
    std::vector<int> cornernodesinallsurfaces = getfacesdefinitionsbasedonnodes();
    std::vector<int> edgesinallsurfaces = getfacesdefinitionsbasedonedges();
    
    if (faceistriangle)
    {
        numberofnodesinsurface = ( (order+1)*(order+1) + (order+1) )/2;
        numberofcornernodesinsurface = 3;
        numberofedgesinsurface = 3;
        offset = 0;
    }
    if (faceisquadrangle)
    {
        numberofnodesinsurface = (order+1)*(order+1);
        numberofcornernodesinsurface = 4;
        numberofedgesinsurface = 4;
        offset = 3*counttriangularfaces();
    }
    int numberofinteriorsurfacenodes = numberofnodesinsurface - numberofedgesinsurface * numberofinteriorlinenodes - numberofcornernodesinsurface;
    std::vector<int> nodesinsurface(numberofnodesinsurface);
    
    // Add the corner nodes of the surface to 'nodesinsurface':
    for (int i = 0; i < numberofcornernodesinsurface; i++)
        nodesinsurface[i] = curvednodelist[cornernodesinallsurfaces[offset+numberofcornernodesinsurface*surfaceindex+i]];

    // Add to 'nodesinsurface' the nodes interior to the edges of the surface (in the correct edge direction):
    int currentedgenumber;
    int index = numberofcornernodesinsurface;
    for (int j = 0; j < numberofedgesinsurface; j++)
    {
        currentedgenumber = edgesinallsurfaces[offset+numberofedgesinsurface*surfaceindex+j];
        // If we follow the edge in the positive defined direction the nodes appear in the correct order in 'curvednodelist':
        if (currentedgenumber > 0)
        {
            // We added 1 to the edge number to be able to give it an orientation sign:
            currentedgenumber = currentedgenumber - 1;
            // Add the interior nodes of the jth edge:
            for (int i = countnodes()+currentedgenumber*numberofinteriorlinenodes; i < countnodes()+(currentedgenumber+1)*numberofinteriorlinenodes; i++)
            {
                nodesinsurface[index] = curvednodelist[i];
                index = index + 1;
            }
        }
        // If we follow the edge in the opposite direction the nodes appear in the reverse order in 'curvednodelist':
        else
        {
            currentedgenumber = -currentedgenumber - 1;
            for (int i = countnodes()+(currentedgenumber+1)*numberofinteriorlinenodes-1; i >= countnodes()+currentedgenumber*numberofinteriorlinenodes; i--)
            {
                nodesinsurface[index] = curvednodelist[i];
                index = index + 1;
            }
        }
    }
    // Add the inner surface nodes of the surface to 'nodesinsurface':

    // In case of a quadrangle surface we have to skip all the interior nodes of the triangular surfaces (if any)
    // since the triangles appear before the quadrangles in the inner face node list of the element object:
    int numberofnodesintriangle = ( (order+1)*(order+1) + (order+1) )/2;
    int numberofinteriortrianglenodes = numberofnodesintriangle - 3 * numberofinteriorlinenodes - 3;
    int toskip;
    if (faceistriangle)
        toskip = 0;
    if (faceisquadrangle)    
        toskip = numberofinteriortrianglenodes*counttriangularfaces();
    
    for (int i=countnodes()+countedges()*numberofinteriorlinenodes+toskip+surfaceindex*numberofinteriorsurfacenodes; i < countnodes()+countedges()*numberofinteriorlinenodes+toskip+(surfaceindex+1)*numberofinteriorsurfacenodes; i++)
    {
        nodesinsurface[index] = curvednodelist[i];
        index = index + 1;
    }

    return nodesinsurface;
}

std::vector<int> element::getedgesdefinitionsbasedonnodes(void)
{
    switch (gettypenumber())
    {
        // Point:
        case 0:
            return {};
        // Line:
        case 1:
            return {0,1};
        // Triangle:
        case 2:
            return {0,1,1,2,2,0};
        // Quadrangle:
        case 3:
            return {0,1,1,2,2,3,3,0};
        // Tetrahedron:
        case 4:
            return {0,1,1,2,2,0,3,0,3,2,3,1};
        // Hexahedron:
        case 5:
            return {0,1,0,3,0,4,1,2,1,5,2,3,2,6,3,7,4,5,4,7,5,6,6,7};
        // Prism:
        case 6:
            return {0,1,0,2,0,3,1,2,1,4,2,5,3,4,3,5,4,5};
        // Pyramid:
        case 7:
            return {0,1,0,3,0,4,1,2,1,4,2,3,2,4,3,4};
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<int> element::getfacesdefinitionsbasedonnodes(void)
{
    switch (gettypenumber())
    {
        // Point:
        case 0:
            return {};
        // Line:
        case 1:
            return {};
        // Triangle:
        case 2:
            return {0,1,2};
        // Quadrangle:
        case 3:
            return {0,1,2,3};
        // Tetrahedron:
        case 4:
            return {0,2,1,0,1,3,0,3,2,3,1,2};
        // Hexahedron:
        case 5:
            return {0,3,2,1,0,1,5,4,0,4,7,3,1,2,6,5,2,3,7,6,4,5,6,7};
        // Prism:
        case 6:
            return {0,2,1,3,4,5,0,1,4,3,0,3,5,2,1,2,5,4};
        // Pyramid:
        case 7:
            return {0,1,4,3,0,4,1,2,4,2,3,4,0,3,2,1};
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<int> element::getfacesdefinitionsbasedonedges(void)
{
    switch (gettypenumber())
    {
        // Point:
        case 0:
            return {};
        // Line:
        case 1:
            return {};
        // Triangle:
        case 2:
            return {1,2,3};
        // Quadrangle:
        case 3:
            return {1,2,3,4};
        // Tetrahedron:
        case 4:
            return {-3,-2,-1,1,-6,4,-4,5,3,6,2,-5};
        // Hexahedron:
        case 5:
            return {2,-6,-4,-1,1,5,-9,-3,3,10,-8,-2,4,7,-11,-5,6,8,-12,-7,9,11,12,-10};
        // Prism:
        case 6:
            return {2,-4,-1,7,9,-8,1,5,-7,-3,3,8,-6,-2,4,6,-9,-5};
        //Pyramid:
        case 7:
            return {1,5,-3,-2,3,-8,4,7,-5,6,8,-7,2,-6,-4,-1};
    }
    
    throw std::runtime_error(""); // fix return warning
}

bool element::iselementedgeorface(void)
{
    int straighttypenumber = gettypenumber();
    // Lines have number 1, triangles 2 and quadrangles 3:
    return (straighttypenumber == 1 || straighttypenumber == 2 || straighttypenumber == 3);
}

std::vector<double> element::listnodecoordinates(void)
{
    std::vector<double> output(3*countcurvednodes(), 0.0);

    int numnodesinline = getcurvatureorder() + 1;
 
     int index = 0;

    switch (gettypenumber())
    {
        // Point:
        case 0:
            output = {0.0,0.0,0.0};
            break;
        // Line:
        case 1:
            for (int i = 0; i < numnodesinline; i++)
                output[3*i+0] = -1.0+2.0/(numnodesinline-1)*i;
            break;
        // Triangle:
        case 2:
            for (int i = 0; i < numnodesinline; i++)
            {
                for (int j = 0; j < numnodesinline-i; j++)
                {
                    output[3*index+0] = 1.0/(numnodesinline-1)*j;
                    output[3*index+1] = 1.0/(numnodesinline-1)*i;
                    index++;
                }
            }
            break;
        // Quadrangle:
        case 3:
            for (int i = 0; i < numnodesinline; i++)
            {
                for (int j = 0; j < numnodesinline; j++)
                {
                    output[3*index+0] = -1.0+2.0/(numnodesinline-1)*j;
                    output[3*index+1] = -1.0+2.0/(numnodesinline-1)*i;
                    index++;
                }
            }
            break;
        // Tetrahedron:
        case 4:
            for (int i = 0; i < numnodesinline; i++)
            {
                for (int j = 0; j < numnodesinline-i; j++)
                {
                    for (int k = 0; k < numnodesinline-i-j; k++)
                    {
                        output[3*index+0] = 1.0/(numnodesinline-1)*k;
                        output[3*index+1] = 1.0/(numnodesinline-1)*j;
                        output[3*index+2] = 1.0/(numnodesinline-1)*i;
                        index++;
                    }
                }
            }
            break;
        // Hexahedron:
        case 5:
            for (int i = 0; i < numnodesinline; i++)
            {
                for (int j = 0; j < numnodesinline; j++)
                {
                    for (int k = 0; k < numnodesinline; k++)
                    {
                        output[3*index+0] = -1.0+2.0/(numnodesinline-1)*k;
                        output[3*index+1] = -1.0+2.0/(numnodesinline-1)*j;
                        output[3*index+2] = -1.0+2.0/(numnodesinline-1)*i;
                        index++;
                    }
                }
            }
            break;
        // Prism:
        case 6:
            for (int i = 0; i < numnodesinline; i++)
            {
                for (int j = 0; j < numnodesinline; j++)
                {
                    for (int k = 0; k < numnodesinline-j; k++)
                    {
                        output[3*index+0] = 1.0/(numnodesinline-1)*k;
                        output[3*index+1] = 1.0/(numnodesinline-1)*j;
                        output[3*index+2] = -1.0+2.0/(numnodesinline-1)*i;
                        index++;
                    }
                }
            }
            break;
        //Pyramid:
        case 7:
            output = {-1.0, -1.0, 0, 1.0, -1.0, 0, -1.0, 1.0, 0.0, 1.0, 1.0, 0, 0.0, 0.0, 1.0};
            if (getcurvatureorder() > 1)
            {             
                logs log;
                log.msg() << "Error in 'element' object: coordinates of order 2 and above not defined for pyramids" << std::endl;
                log.error();
            }
            break;
    }

    return output;
}

int element::deducetypenumber(int elemdim, int numnodes)
{
    switch (numnodes)
    {
        // Point:
        case 1:
            return 0;
        // Line:
        case 2:
            return 1;
        // Triangle:
        case 3:
            return 2;
        // Quadrangle or tetrahedron:
        case 4:
        {
            if (elemdim == 2)
                return 3;
            if (elemdim == 3)
                return 4;
            break;
        }
        // Pyramid:
        case 5:
            return 7;
        // Prism:
        case 6:
            return 6;
        // Hexahedron:
        case 8:
            return 5;
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<double> element::calculatecoordinates(std::vector<double>& refcoords, std::vector<double>& nodecoords, int fi, bool returnnodecoords)
{
    // Need to be ultrafast for straight elements in mesh adaptivity
    if (returnnodecoords == false)
    {
        if (curvedtypenumber == 1)
        {
            int numrefs = refcoords.size()/3;
            std::vector<double> outfast(3*numrefs);
            for (int i = 0; i < numrefs; i++)
            {
                double ki = refcoords[3*i+0];
                
                double ff0 = +0.5-0.5*ki;
                double ff1 = +0.5+0.5*ki;
                
                outfast[3*i+0] = nodecoords[fi+3*0+0]*ff0 + nodecoords[fi+3*1+0]*ff1;
                outfast[3*i+1] = nodecoords[fi+3*0+1]*ff0 + nodecoords[fi+3*1+1]*ff1;
                outfast[3*i+2] = nodecoords[fi+3*0+2]*ff0 + nodecoords[fi+3*1+2]*ff1;
            }
            return outfast;
        }
        
        if (curvedtypenumber == 2)
        {
            int numrefs = refcoords.size()/3;
            std::vector<double> outfast(3*numrefs);
            for (int i = 0; i < numrefs; i++)
            {
                double ki = refcoords[3*i+0];
                double eta = refcoords[3*i+1];
                
                double ff0 = +1-1*eta-1*ki;
                double ff1 = +1*ki;
                double ff2 = +1*eta;
                
                outfast[3*i+0] = nodecoords[fi+3*0+0]*ff0 + nodecoords[fi+3*1+0]*ff1 + nodecoords[fi+3*2+0]*ff2;
                outfast[3*i+1] = nodecoords[fi+3*0+1]*ff0 + nodecoords[fi+3*1+1]*ff1 + nodecoords[fi+3*2+1]*ff2;
                outfast[3*i+2] = nodecoords[fi+3*0+2]*ff0 + nodecoords[fi+3*1+2]*ff1 + nodecoords[fi+3*2+2]*ff2;
            }
            return outfast;
        }

        if (curvedtypenumber == 3)
        {
            int numrefs = refcoords.size()/3;
            std::vector<double> outfast(3*numrefs);
            for (int i = 0; i < numrefs; i++)
            {
                double ki = refcoords[3*i+0];
                double eta = refcoords[3*i+1];
                
                double ff0 = +0.25-0.25*eta-0.25*ki+0.25*ki*eta;
                double ff1 = +0.25-0.25*eta+0.25*ki-0.25*ki*eta;
                double ff2 = +0.25+0.25*eta+0.25*ki+0.25*ki*eta;
                double ff3 = +0.25+0.25*eta-0.25*ki-0.25*ki*eta;

                outfast[3*i+0] = nodecoords[fi+3*0+0]*ff0 + nodecoords[fi+3*1+0]*ff1 + nodecoords[fi+3*2+0]*ff2 + nodecoords[fi+3*3+0]*ff3;
                outfast[3*i+1] = nodecoords[fi+3*0+1]*ff0 + nodecoords[fi+3*1+1]*ff1 + nodecoords[fi+3*2+1]*ff2 + nodecoords[fi+3*3+1]*ff3;
                outfast[3*i+2] = nodecoords[fi+3*0+2]*ff0 + nodecoords[fi+3*1+2]*ff1 + nodecoords[fi+3*2+2]*ff2 + nodecoords[fi+3*3+2]*ff3;
            }
            return outfast;
        }
        
        if (curvedtypenumber == 4)
        {
            int numrefs = refcoords.size()/3;
            std::vector<double> outfast(3*numrefs);
            for (int i = 0; i < numrefs; i++)
            {
                double ki = refcoords[3*i+0];
                double eta = refcoords[3*i+1];
                double phi = refcoords[3*i+2];
                
                double ff0 = +1-1*phi-1*eta-1*ki;
                double ff1 = +1*ki;
                double ff2 = +1*eta;
                double ff3 = +1*phi;
                
                outfast[3*i+0] = nodecoords[fi+3*0+0]*ff0 + nodecoords[fi+3*1+0]*ff1 + nodecoords[fi+3*2+0]*ff2 + nodecoords[fi+3*3+0]*ff3;
                outfast[3*i+1] = nodecoords[fi+3*0+1]*ff0 + nodecoords[fi+3*1+1]*ff1 + nodecoords[fi+3*2+1]*ff2 + nodecoords[fi+3*3+1]*ff3;
                outfast[3*i+2] = nodecoords[fi+3*0+2]*ff0 + nodecoords[fi+3*1+2]*ff1 + nodecoords[fi+3*2+2]*ff2 + nodecoords[fi+3*3+2]*ff3;
            }
            return outfast;
        }
    }


    int numpolys = mypolynomials.count();
    if (numpolys > 0)
    {
        if (returnnodecoords)
        {
            std::vector<double> out(3*numpolys);
            for (int i = 0; i < 3*numpolys; i++)
                out[i] = nodecoords[fi+i];
            return out;
        }
    
        int numrefs = refcoords.size()/3;
        
        std::vector<double> sf, evaluationpoint;
        std::vector<double> output(3*numrefs, 0.0);        

        for (int i = 0; i < numrefs; i++)
        {
            evaluationpoint = {refcoords[3*i+0],refcoords[3*i+1],refcoords[3*i+2]};
            
            mypolynomials.evalatsingle(evaluationpoint, sf);

            for (int c = 0; c < numpolys; c++)
            {
                output[3*i+0] += nodecoords[fi + 3*c+0] * sf[c];
                output[3*i+1] += nodecoords[fi + 3*c+1] * sf[c];
                output[3*i+2] += nodecoords[fi + 3*c+2] * sf[c];
            }
        }
        return output;
    }
    else
    {
        lagrangeformfunction lff(gettypenumber(), getcurvatureorder(), {});
        mypolynomials = polynomials(lff.getformfunctionpolynomials());
        return calculatecoordinates(refcoords, nodecoords, fi, returnnodecoords);
    }
}

std::vector<int> element::fullsplitcount(int n)
{
    std::vector<int> numsons = {1,2,4,4,8,8,8,6};
    std::vector<int> output(8,0);

    int typenum = gettypenumber();
    
    output[typenum] = std::pow(numsons[typenum],n);
    
    // Pyramids also give tetrahedra when split:
    if (typenum == 7)
        output[4] = std::pow(2,n+1)*(std::pow(4,n)-std::pow(3,n));
    
    return output;
}

void element::fullsplit(int n, std::vector<std::vector<double>>& splitcoords, std::vector<double>& unsplitcoords)
{
    if (n == 0)
    {
        splitcoords = std::vector<std::vector<double>>(8,std::vector<double>(0));
        splitcoords[gettypenumber()] = unsplitcoords;
        return;
    }
    // Recursive call:
    if (n > 1)
    {
        std::vector<std::vector<double>> cursplitcoords;
        fullsplit(1, cursplitcoords, unsplitcoords);
        fullsplit(n-1, splitcoords, cursplitcoords[gettypenumber()]);
        // Treat the tetrahedra from the split pyramids:
        if (gettypenumber() == 7)
        {
            element mytet(4,getcurvatureorder());
            std::vector<std::vector<double>> tetsplitcoords;
            mytet.fullsplit(n-1, tetsplitcoords, cursplitcoords[4]);
            int cursize = splitcoords[4].size();
            int sizetoadd = tetsplitcoords[4].size();
            splitcoords[4].resize(cursize+sizetoadd);
            for (int i = 0; i < sizetoadd; i++)
                splitcoords[4][cursize+i] = tetsplitcoords[4][i];
        }
        return;
    }
    
    int tn = gettypenumber();
    int co = getcurvatureorder();
    int nn = countnodes();
    int ncn = countcurvednodes();
    int ne = unsplitcoords.size()/3/ncn;
    std::vector<int> splitcount = fullsplitcount(1);
    int ns = splitcount[tn];
    element straighelem(tn);
    
    lagrangeformfunction lff(tn, co, {});
    std::vector<double> curvedrefcoords = lff.getnodecoordinates();
    
    // Preallocate:
    splitcoords = std::vector<std::vector<double>>(8,std::vector<double>(0));
    splitcoords[tn] = std::vector<double>(3*ncn*ns*ne);
    // Define only once:
    std::vector< std::vector<std::vector<double>> > cornerrefcoords(3);
    int numcases = 1;
    if (tn == 4)
        numcases = 3;
    for (int c = 0; c < numcases; c++)
        fullsplit(cornerrefcoords[c], c);
        
    std::vector< std::vector<std::vector<double>> > curvedsubcoords(numcases, std::vector<std::vector<double>>(ns));
    for (int i = 0; i < ns; i++)
    {
        for (int c = 0; c < numcases; c++)
        {
            std::vector<double> cc(3*nn);
            for (int j = 0; j < 3*nn; j++)
                cc[j] = cornerrefcoords[c][tn][3*nn*i+j];
            curvedsubcoords[c][i] = straighelem.calculatecoordinates(curvedrefcoords, cc);
        }
    }
    
    // Populate:
    int index = 0;
    for (int e = 0; e < ne; e++)
    {
        std::vector<double> coords(3*ncn);
        for (int i = 0; i < 3*ncn; i++)
            coords[i] = unsplitcoords[3*ncn*e+i];
    
        int throughedgenum = 0;
        if (tn == 4)
            throughedgenum = choosethroughedge(coords);
        
        for (int i = 0; i < ns; i++)
        {
            std::vector<double> cursplit = calculatecoordinates(curvedsubcoords[throughedgenum][i], coords);

            for (int j = 0; j < cursplit.size(); j++)
                splitcoords[tn][index+j] = cursplit[j];
            index += cursplit.size();
        }
    }
    
    // Also add the tetrahedra created during a pyramid split:
    if (tn == 7)
    {
        element mystraighttet(4);
        lagrangeformfunction lfftet(4, co, {});
        std::vector<double> curvedtetrefcoords = lfftet.getnodecoordinates();
    
        splitcoords[4] = std::vector<double>(curvedtetrefcoords.size()*4*ne);
    
        std::vector<std::vector<double>> curvedtetsubcoords(4);
        for (int i = 0; i < 4; i++)
        {
            std::vector<double> cc(12);
            for (int j = 0; j < 12; j++)
                cc[j] = cornerrefcoords[0][4][12*i+j];
            curvedtetsubcoords[i] = mystraighttet.calculatecoordinates(curvedtetrefcoords, cc);
        }
    
        int index = 0;
        for (int e = 0; e < ne; e++)
        {
            std::vector<double> coords(3*ncn);
            for (int i = 0; i < 3*ncn; i++)
                coords[i] = unsplitcoords[3*ncn*e+i];
        
            for (int i = 0; i < 4; i++)
            {
                std::vector<double> cursplit = calculatecoordinates(curvedtetsubcoords[i], coords);

                for (int j = 0; j < cursplit.size(); j++)
                    splitcoords[4][index+j] = cursplit[j];
                index += cursplit.size();
            }
        }
    }
    
}

void element::fullsplit(std::vector<std::vector<double>>& cornerrefcoords, int throughedgenum)
{
    cornerrefcoords = std::vector<std::vector<double>>(8,std::vector<double>(0));

    switch (gettypenumber())
    {
        case 0:
        {
            cornerrefcoords[0] = {0,0,0};
            break;
        }
        case 1:
        {
            cornerrefcoords[1] = {-1,0,0,0,0,0, 0,0,0,1,0,0};
            break;
        }
        case 2:
        {
            cornerrefcoords[2] = {0,0,0,0.5,0,0,0,0.5,0, 0.5,0,0,0.5,0.5,0,0,0.5,0, 0.5,0,0,1,0,0,0.5,0.5,0, 0,0.5,0,0.5,0.5,0,0,1,0};
            break;
        }
        case 3:
        {
            cornerrefcoords[3] = {-1,-1,0,0,-1,0,0,0,0,-1,0,0, 0,-1,0,1,-1,0,1,0,0,0,0,0, -1,0,0,0,0,0,0,1,0,-1,1,0, 0,0,0,1,0,0,1,1,0,0,1,0};
            break;
        }
        case 4:
        {
            if (throughedgenum == 0)
                cornerrefcoords[4] = {0.5,0,0,0,0.5,0.5,0,0,0.5,0.5,0,0.5, 0.5,0,0,0,0.5,0,0,0,0.5,0,0.5,0.5, 0.5,0,0,0.5,0.5,0,0,0.5,0.5,0.5,0,0.5, 0.5,0,0,0.5,0.5,0,0,0.5,0,0,0.5,0.5, 0,0,1,0,0,0.5,0,0.5,0.5,0.5,0,0.5, 0,1,0,0,0.5,0,0.5,0.5,0,0,0.5,0.5, 1,0,0,0.5,0.5,0,0.5,0,0,0.5,0,0.5, 0,0,0,0.5,0,0,0,0.5,0,0,0,0.5};
            if (throughedgenum == 1)
                cornerrefcoords[4] = {0.5,0.5,0,0,0.5,0.5,0,0,0.5,0.5,0,0.5, 0.5,0.5,0,0,0.5,0,0,0,0.5,0,0.5,0.5, 0.5,0,0,0.5,0.5,0,0,0,0.5,0.5,0,0.5, 0.5,0,0,0.5,0.5,0,0,0.5,0,0,0,0.5, 0,0,1,0,0,0.5,0,0.5,0.5,0.5,0,0.5, 0,1,0,0,0.5,0,0.5,0.5,0,0,0.5,0.5, 1,0,0,0.5,0.5,0,0.5,0,0,0.5,0,0.5, 0,0,0,0.5,0,0,0,0.5,0,0,0,0.5};
            if (throughedgenum == 2)
                cornerrefcoords[4] = {0,0.5,0,0,0.5,0.5,0,0,0.5,0.5,0,0.5, 0.5,0.5,0,0,0.5,0.5,0,0.5,0,0.5,0,0.5, 0.5,0,0,0,0.5,0,0,0,0.5,0.5,0,0.5, 0.5,0,0,0.5,0.5,0,0,0.5,0,0.5,0,0.5, 0,0,1,0,0,0.5,0,0.5,0.5,0.5,0,0.5, 0,1,0,0,0.5,0,0.5,0.5,0,0,0.5,0.5, 1,0,0,0.5,0.5,0,0.5,0,0,0.5,0,0.5, 0,0,0,0.5,0,0,0,0.5,0,0,0,0.5};
            break;
        }
        case 5:
        {
            cornerrefcoords[5] = {-1,-1,-1,0,-1,-1,0,0,-1,-1,0,-1,-1,-1,0,0,-1,0,0,0,0,-1,0,0, 0,-1,-1,1,-1,-1,1,0,-1,0,0,-1,0,-1,0,1,-1,0,1,0,0,0,0,0, -1,0,-1,0,0,-1,0,1,-1,-1,1,-1,-1,0,0,0,0,0,0,1,0,-1,1,0, 0,0,-1,1,0,-1,1,1,-1,0,1,-1,0,0,0,1,0,0,1,1,0,0,1,0, -1,-1,0,0,-1,0,0,0,0,-1,0,0,-1,-1,1,0,-1,1,0,0,1,-1,0,1, 0,-1,0,1,-1,0,1,0,0,0,0,0,0,-1,1,1,-1,1,1,0,1,0,0,1, -1,0,0,0,0,0,0,1,0,-1,1,0,-1,0,1,0,0,1,0,1,1,-1,1,1, 0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1};
            break;
        }
        case 6:
        {
            cornerrefcoords[6] = {0,0,-1,0.5,0,-1,0,0.5,-1,0,0,0,0.5,0,0,0,0.5,0, 0.5,0,-1,0.5,0.5,-1,0,0.5,-1,0.5,0,0,0.5,0.5,0,0,0.5,0, 0.5,0,-1,1,0,-1,0.5,0.5,-1,0.5,0,0,1,0,0,0.5,0.5,0, 0,0.5,-1,0.5,0.5,-1,0,1,-1,0,0.5,0,0.5,0.5,0,0,1,0, 0,0,0,0.5,0,0,0,0.5,0,0,0,1,0.5,0,1,0,0.5,1, 0.5,0,0,0.5,0.5,0,0,0.5,0,0.5,0,1,0.5,0.5,1,0,0.5,1, 0.5,0,0,1,0,0,0.5,0.5,0,0.5,0,1,1,0,1,0.5,0.5,1, 0,0.5,0,0.5,0.5,0,0,1,0,0,0.5,1,0.5,0.5,1,0,1,1};
            break;
        }
        case 7:
        {
            cornerrefcoords[4] = {0,-1,0,-1.0,-1.0,0.5,1.0,-1.0,0.5,0,0,0, 1.0,-1.0,0.5,1.0,1.0,0.5,1,0,0,0,0,0, 1.0,1.0,0.5,-1.0,1.0,0.5,0,1,0,0,0,0, -1,0,0,-1.0,1.0,0.5,-1.0,-1.0,0.5,0,0,0};
            cornerrefcoords[7] = {-1.0,-1.0,0.5,1.0,-1.0,0.5,1.0,1.0,0.5,-1.0,1.0,0.5,0,0,1, -1.0,-1.0,0.5,-1.0,1.0,0.5,1.0,1.0,0.5,1.0,-1.0,0.5,0,0,0, -1,-1,0,0,-1,0,0,0,0,-1,0,0,-1.0,-1.0,0.5, -1,0,0,0,0,0,0,1,0,-1,1,0,-1.0,1.0,0.5, 0,-1,0,1,-1,0,1,0,0,0,0,0,1.0,-1.0,0.5, 0,0,0,1,0,0,1,1,0,0,1,0,1.0,1.0,0.5};
            break;
        }
    }
}

int element::choosethroughedge(std::vector<double>& nodecoords)
{
    // Mid-edge reference coordinates:
    std::vector<double> refcoords = {0.5,0,0, 0.5,0.5,0, 0,0.5,0, 0,0,0.5, 0,0.5,0.5, 0.5,0,0.5};

    std::vector<double> calced = calculatecoordinates(refcoords, nodecoords);
    
    double l04 = geotools::getdistance(0,4, calced);
    double l13 = geotools::getdistance(1,3, calced);
    double l25 = geotools::getdistance(2,5, calced);

    // Select shortest edge:
    if (l04 <= l13 && l04 <= l25)
        return 0;
    if (l13 <= l04 && l13 <= l25)
        return 1;
    if (l25 <= l04 && l25 <= l13)
        return 2;
        
    throw std::runtime_error(""); // fix return warning
}

std::vector<std::vector<int>> element::split(int splitnum, std::vector<int>& edgenumbers)
{
    if (gettypenumber() > 4)
    {
        logs log;
        log.msg() << "Error in 'element' object: transition splits not defined yet for " << gettypenameconjugation(2) << std::endl;
        log.error();
    }

    switch (gettypenumber())
    {
        case 0:
            return {{0},{},{},{},{},{},{},{}};
        case 1:
            return splitline(splitnum);
        case 2:
            return splittriangle(splitnum, edgenumbers);
        case 3:
            return splitquadrangle(splitnum);
        case 4:
        {
            // Too slow to be called each time. Compute and store.
            std::vector<std::vector<int>> splitdef;
            if (universe::getsplitdefinition(splitdef, 4, splitnum, edgenumbers))
                return splitdef;
                
            splitdef = splittetrahedron(splitnum, edgenumbers);
            universe::setsplitdefinition(splitdef, 4, splitnum, edgenumbers);
            return splitdef;
        }
    }
    
    throw std::runtime_error(""); // fix return warning
}
        
std::vector<std::vector<int>> element::splitline(int splitnum)
{
    if (splitnum == 0)
        return {{},{0,1},{},{},{},{},{},{}};
    else
        return {{},{0,2, 2,1},{},{},{},{},{},{}};
}
        
std::vector<std::vector<int>> element::splittriangle(int splitnum, std::vector<int>& edgenumbers)
{
    switch (splitnum)
    {
        case 0:
            return {{},{},{0,1,2},{},{},{},{},{}};
        case 1:
            return {{},{},{0,1,5, 5,1,2},{},{},{},{},{}};
        case 2:
            return {{},{},{0,1,4, 0,4,2},{},{},{},{},{}};
        case 3:
        {
            if (edgenumbers[2] < edgenumbers[1])
                return {{},{},{0,1,5, 5,1,4, 5,4,2},{},{},{},{},{}};
            else
                return {{},{},{0,1,4, 0,4,5, 5,4,2},{},{},{},{},{}};
        }
        case 4:
            return {{},{},{0,3,2, 3,1,2},{},{},{},{},{}};
        case 5:
        {
            if (edgenumbers[2] < edgenumbers[0])
                return {{},{},{0,3,5, 5,3,1, 5,1,2},{},{},{},{},{}};
            else
                return {{},{},{0,3,5, 5,3,2, 3,1,2},{},{},{},{},{}};
        }
        case 6:
        {
            if (edgenumbers[1] < edgenumbers[0])
                return {{},{},{0,3,4, 0,4,2, 3,1,4},{},{},{},{},{}};
            else
                return {{},{},{0,3,2, 3,4,2, 3,1,4},{},{},{},{},{}};
        }
        case 7:
            return {{},{},{0,3,5, 3,1,4, 3,4,5, 5,4,2},{},{},{},{},{}};
    }
    
    throw std::runtime_error(""); // fix return warning
}    

std::vector<std::vector<int>> element::splitquadrangle(int splitnum)
{
    switch (splitnum)
    {
        case 0:
            return {{},{},{},{0,1,2,3},{},{},{},{}};
        case 1:
            return {{},{},{0,1,7, 7,1,2, 7,2,3},{},{},{},{},{}};
        case 2:
            return {{},{},{0,6,3, 0,1,6, 1,2,6},{},{},{},{},{}};
        case 3:
            return {{},{},{},{0,1,8,7, 7,8,6,3, 8,1,2,6},{},{},{},{}};
        case 4:
            return {{},{},{0,1,5, 0,5,3, 3,5,2},{},{},{},{},{}};
        case 5:
            return {{},{},{},{0,1,5,7, 7,5,2,3},{},{},{},{}};
        case 6:
            return {{},{},{},{0,1,5,8, 0,8,6,3, 8,5,2,6},{},{},{},{}};
        case 7:
            return {{},{},{7,5,6, 7,6,3, 5,2,6},{0,1,5,7},{},{},{},{}};
        case 8:
            return {{},{},{0,4,3, 4,2,3, 4,1,2},{},{},{},{},{}};
        case 9:
            return {{},{},{},{0,4,8,7, 4,1,2,8, 7,8,2,3},{},{},{},{}};
        case 10:
            return {{},{},{},{0,4,6,3, 4,1,2,6},{},{},{},{}};
        case 11:
            return {{},{},{0,4,7, 4,6,7, 7,6,3},{4,1,2,6},{},{},{},{}};
        case 12:
            return {{},{},{},{0,4,8,3, 4,1,5,8, 8,5,2,3},{},{},{},{}};
        case 13:
            return {{},{},{0,4,7, 7,4,5, 4,1,5},{7,5,2,3},{},{},{},{}};
        case 14:
            return {{},{},{4,1,5, 4,5,6, 5,2,6},{0,4,6,3},{},{},{},{}};
        case 15:
            return {{},{},{},{0,4,8,7, 4,1,5,8, 7,8,6,3, 8,5,2,6},{},{},{},{}};
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<std::vector<int>> element::splittetrahedron(int splitnum, std::vector<int>& edgenumbers)
{
    std::vector<bool> isedgecut = gentools::inttobinary(6, splitnum);

    // The triangle element object can be straight or curved (no impact here):
    element mytri(2);

    // Split the faces and compute their node connectivity:
    std::vector<std::vector<bool>> faceconnectivity(4);
    std::vector<int> faceedgedef = {2,1,0, 0,5,3, 3,4,2, 5,1,4};
    for (int i = 0; i < 4; i++)
    {
        int ea = faceedgedef[3*i+0];
        int eb = faceedgedef[3*i+1];
        int ec = faceedgedef[3*i+2];

        std::vector<bool> iec = {isedgecut[ea], isedgecut[eb], isedgecut[ec]};
        std::vector<int> triedgenums = {edgenumbers[ea], edgenumbers[eb], edgenumbers[ec]};
        std::vector<std::vector<int>> trisplitdefinition = mytri.split(gentools::binarytoint(iec), triedgenums);
        mytri.getsplitconnectivity(faceconnectivity[i], trisplitdefinition);
    }
    // Compute the tetrahedron node connectivity:
    std::vector<bool> connectivity;
    getsplitconnectivity(connectivity, faceconnectivity);

    
    // Total number of edges cut:
    int countcuts = 0;
    for (int i = 0; i < 6; i++)
    {
        if (isedgecut[i])
            countcuts++;
    }

    // In case a through-edge has to be inserted:
    std::vector<int> facingedges = {0,4, 1,3, 2,5};

    std::vector<int> edgefacingpair = {};
    for (int i = 0; i < 3; i++)
    {
        if (isedgecut[facingedges[2*i+0]] && isedgecut[facingedges[2*i+1]])
            edgefacingpair.push_back(i);
    }
    
    std::vector<std::vector<bool>> connectivitythroughedge(edgefacingpair.size());
    
    // Add connection:
    for (int i = 0; i < edgefacingpair.size(); i++)
    {
        connectivitythroughedge[i] = connectivity;
    
        int m = 4 + facingedges[2*edgefacingpair[i]+0];
        int n = 4 + facingedges[2*edgefacingpair[i]+1];

        connectivitythroughedge[i][m*10+n] = true;
        connectivitythroughedge[i][n*10+m] = true;
    }

    // Connect the tets:
    std::vector<std::vector<bool>> tetdefs;
    deducetets(connectivity, tetdefs);
    int numtets = tetdefs.size();
    
    std::vector<int> minnumtet = {1,2,3,4,5,7,8};
    // If a through-edge must be inserted:
    int ind = 0;
    while (numtets < minnumtet[countcuts])
    {
        deducetets(connectivitythroughedge[ind], tetdefs);
        numtets = tetdefs.size();
        ind++;
    }

    // Convert the boolean tetrahedra definition to node numbers:
    std::vector<std::vector<int>> output(8, std::vector<int>(0));
    output[4] = std::vector<int>(numtets*4);
    for (int i = 0; i < numtets; i++)
    {
        int index = 0;
        for (int j = 0; j < 10; j++)
        {
            if (tetdefs[i][j])
            {
                output[4][4*i+index] = j;
                index++;   
            }
        }
    }
    
    // Reorder the tetrahedra nodes if needed:
    std::vector<double> rc = {0,0,0, 1,0,0, 0,1,0, 0,0,1, 0.5,0,0, 0.5,0.5,0, 0,0.5,0, 0,0,0.5, 0,0.5,0.5, 0.5,0,0.5};
    for (int i = 0; i < numtets; i++)
    {
        int p0 = output[4][4*i+0];
        int p1 = output[4][4*i+1];
        int p2 = output[4][4*i+2];
        int p3 = output[4][4*i+3];
        
        // Perpendicular to bottom tet face:
        std::vector<double> a = {rc[3*p1+0]-rc[3*p0+0], rc[3*p1+1]-rc[3*p0+1], rc[3*p1+2]-rc[3*p0+2]};
        std::vector<double> b = {rc[3*p2+0]-rc[3*p0+0], rc[3*p2+1]-rc[3*p0+1], rc[3*p2+2]-rc[3*p0+2]};
        std::vector<double> c = {rc[3*p3+0]-rc[3*p0+0], rc[3*p3+1]-rc[3*p0+1], rc[3*p3+2]-rc[3*p0+2]};
        // Cross product:
        std::vector<double> cp = {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
        // Check if p3 is above or below (if below the surface orientation is flipped):
        bool isbelow = true;
        // Cannot be zero for non-degenerate tet:
        if (cp[0]*c[0]+cp[1]*c[1]+cp[2]*c[2] > 0)
            isbelow = false;
            
        if (isbelow)
        {
            output[4][4*i+1] = p2;
            output[4][4*i+2] = p1;
        }
    }
    
    return output;
}

void element::getsplitconnectivity(std::vector<bool>& connectivity, std::vector<std::vector<int>>& splitdefinition)
{
    connectivity = std::vector<bool>(6*6, false);

    for (int j = 0; j < splitdefinition[2].size()/3; j++)
    {
        for (int n = 0; n < 3; n++)
        {
            int curnode = splitdefinition[2][3*j+n];
            int prevnode = splitdefinition[2][3*j+(n+3-1)%3];
            
            connectivity[6*curnode+curnode] = true;
            connectivity[6*prevnode+curnode] = true;
            connectivity[6*curnode+prevnode] = true;
        }
    }
}

void element::getsplitconnectivity(std::vector<bool>& volumeconnectivity, std::vector<std::vector<bool>>& faceconnectivity)
{
    volumeconnectivity = std::vector<bool>(10*10, false);
    
    std::vector<int> fnd = getfacesdefinitionsbasedonnodes();
    std::vector<int> fed = getfacesdefinitionsbasedonedges();
    // Face edge definition is provided in a special format:
    for (int i = 0; i < fed.size(); i++)
        fed[i] = std::abs(fed[i])-1;                                        
    
    // Loop on each face:
    for (int f = 0; f < 4; f++)
    {
        // Edge nodes come after the 4 corner nodes:
        std::vector<int> nmap = {fnd[3*f+0], fnd[3*f+1], fnd[3*f+2], 4+fed[3*f+0], 4+fed[3*f+1], 4+fed[3*f+2]};
    
        for (int m = 0; m < 6; m++)
        {
            for (int n = 0; n < 6; n++)
                volumeconnectivity[10*nmap[m]+nmap[n]] = faceconnectivity[f][6*m+n];
        }
    }
}

void element::deducetets(std::vector<bool>& connectivity, std::vector<std::vector<bool>>& tetdefs, int originnode, int numinloop, int curnode, std::vector<bool> isnodeused)
{
    // Loop on all connected nodes (except the node itself and previously used nodes):
    for (int i = 0; i < 10; i++)
    {
        if (i == curnode || connectivity[10*curnode+i] == false)
            continue;
            
        if (numinloop < 4 && isnodeused[i])
            continue;
            
        // The new node must be connected directly to all other ones:
        bool isconnectedtoall = true;
        for (int j = 0; j < 10; j++)
        {
            if (isnodeused[j] && connectivity[10*i+j] == false)
            {
                isconnectedtoall = false;
                break;
            }
        }
        if (isconnectedtoall == false)
            continue;

        // Terminate:
        if (numinloop == 4)
        {
            // New tet found:
            if (i == originnode)
                tetdefs.push_back(isnodeused);
            continue;
        }
        std::vector<bool> isused = isnodeused;
        isused[i] = true;
        deducetets(connectivity, tetdefs, originnode, numinloop+1, i, isused);
    }
}

void element::deducetets(std::vector<bool>& connectivity, std::vector<std::vector<bool>>& tetdefs)
{
    tetdefs = {};
    
    // Each node is a seed for a tet:
    for (int t = 0; t < 10; t++)
    {
        if (connectivity[10*t+t] == false)
            continue;
    
        std::vector<bool> isnodeused(10, false);
        isnodeused[t] = true;
        deducetets(connectivity, tetdefs, t, 1, t, isnodeused);
    }
    
    // Remove duplicates:
    std::sort( tetdefs.begin(), tetdefs.end() );
    tetdefs.erase( std::unique( tetdefs.begin(), tetdefs.end() ), tetdefs.end() );
}

void element::numstorefcoords(std::vector<int>& nums, std::vector<double>& refcoords)
{
    refcoords = std::vector<double>(3*nums.size());

    if (nums.size() == 0)
        return;
        
    std::vector<double> refs;
    
    switch (gettypenumber())
    {
        case 0:
        {
            refs = {0,0,0};
            break;
        }
        case 1:
        {
            refs = {-1,0,0, 1,0,0, 0,0,0};
            break;
        }
        case 2:
        {
            refs = {0,0,0, 1,0,0, 0,1,0, 0.5,0,0, 0.5,0.5,0, 0,0.5,0};
            break;
        }
        case 3:
        {
            refs = {-1,-1,0, 1,-1,0, 1,1,0, -1,1,0, 0,-1,0, 1,0,0, 0,1,0, -1,0,0, 0,0,0};
            break;
        }
        case 4:
        {
            refs = {0,0,0, 1,0,0, 0,1,0, 0,0,1, 0.5,0,0, 0.5,0.5,0, 0,0.5,0, 0,0,0.5, 0,0.5,0.5, 0.5,0,0.5};
            break;
        }
        case 5:
        {
            refs = {};
            break;
        }
        case 6:
        {
            refs = {};
            break;
        }
        case 7:
        {
            refs = {};
            break;
        }
    }
        
    for (int i = 0; i < nums.size(); i++)
    {
        refcoords[3*i+0] = refs[3*nums[i]+0];
        refcoords[3*i+1] = refs[3*nums[i]+1];
        refcoords[3*i+2] = refs[3*nums[i]+2];
    }
}

void element::write(std::string filename, std::vector<double> coords)
{
    int ncn = countcurvednodes();
    int numelems = coords.size()/3/ncn;
    
    densemat xcoords(numelems, ncn);
    densemat ycoords(numelems, ncn);
    densemat zcoords(numelems, ncn);
    
    densemat vals(numelems, ncn);

    double* xptr = xcoords.getvalues();
    double* yptr = ycoords.getvalues();
    double* zptr = zcoords.getvalues();
    double* vptr = vals.getvalues();
    
    for (int e = 0; e < numelems; e++)
    {
        for (int n = 0; n < ncn; n++)
        {
            xptr[e*ncn+n] = coords[3*ncn*e+3*n+0];
            yptr[e*ncn+n] = coords[3*ncn*e+3*n+1];
            zptr[e*ncn+n] = coords[3*ncn*e+3*n+2];
            vptr[e*ncn+n] = e;
        }
    }
    
    iodata datatowrite(getcurvatureorder(), getcurvatureorder(), true, {});
    
    datatowrite.addcoordinates(gettypenumber(), xcoords, ycoords, zcoords);
    datatowrite.adddata(gettypenumber(), {vals});
    
    iointerface::writetofile(filename, datatowrite);
}

}

