#include "element.h"
#include "geotools.h"
#include "lagrangeformfunction.h"


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
        std::cout << "Error in 'element' object: trying to use undefined element name: " << elementname << std::endl << "Make sure everything is lower case" << std::endl;
        abort();
    }
}

element::element(int number)
{
    if (number < 0)
    {
        std::cout << "Error in 'element' object: cannot define negative element type number " << number << std::endl;
        abort();
    }
    curvedtypenumber = number;
}

element::element(int number, int curvatureorder)
{
    if (number < 0)
    {
        std::cout << "Error in 'element' object: can not define a negative element type number" << std::endl;
        std::cout << "Element type number is " << number << " with " << curvatureorder << " curvature order" << std::endl;
        abort();
    }
    if (curvatureorder <= 0)
    {
        std::cout << "Error in 'element' object: can not define a negative or 0 curvature order" << std::endl;
        std::cout << "Element type number is " << number << " with " << curvatureorder << " curvature order" << std::endl;
        abort();
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
        std::cout << "Error: element type has not been defined yet" << std::endl;
        abort();
    }
    if (nodelist.size() != countcurvednodes())
    {
        std::cout << "Error: trying to define an order " << getcurvatureorder() << " " << gettypename() << " with " << nodelist.size() << " nodes. There should be " << countcurvednodes() << std::endl;
        abort();
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
            curvednumberofnodes = ( pow(order+1,2) + (order+1) )/2;
            break;
        // Quadrangle:
        case 3:
            curvednumberofnodes = pow(order + 1,2);
            break;
        // Tetrahedron:
        case 4:
            curvednumberofnodes = ( (order+1)*(order+2)*(order+3) )/6;
            break;
        // Hexahedron:
        case 5:
            curvednumberofnodes = pow(order + 1,3);
            break;
        // Prism:
        case 6:
            curvednumberofnodes = ( (order + 1)*( pow(order+1,2) + (order+1) ) )/2;
            break;
        // Pyramid:
        case 7:
            // sum of all parallel quadrangles:
            curvednumberofnodes = 0;
            for (int i = 0; i < order + 1; i++)
                curvednumberofnodes = curvednumberofnodes + pow(i + 1,2);
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
    if (curvednodelist.size() != countcurvednodes())
    {
        std::cout << "Error: trying to use an order " << getcurvatureorder() << " " << gettypename() << " with " << curvednodelist.size() << " nodes. There should be " << countcurvednodes() << ". Have you correctly defined the node list?" << std::endl;
        abort();
    }
    if (lineindex+1 > countedges())
    {
        std::cout << "Error: trying to get the " << lineindex+1 << "th line in a " << gettypename() << " but there are only " << countedges() << std::endl;
        abort();
    }
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
    if (curvednodelist.size() != countcurvednodes())
    {
        std::cout << "Error: trying to use an order " << getcurvatureorder() << " " << gettypename() << " with " << curvednodelist.size() << " nodes. There should be " << countcurvednodes() << ". Have you correctly defined the node list?" << std::endl;
        abort();
    }
    if (triangleindex+1 > counttriangularfaces())
    {
        std::cout << "Error: trying to get the " << triangleindex+1 << "th triangle in a " << gettypename() << " but there are only " << counttriangularfaces() << std::endl;
        abort();
    }

    return getnodesinsurface(triangleindex, true, false);
}

std::vector<int> element::getnodesinquadrangle(int quadrangleindex)
{
    if (curvednodelist.size() != countcurvednodes())
    {
        std::cout << "Error: trying to use an order " << getcurvatureorder() << " " << gettypename() << " with " << curvednodelist.size() << " nodes. There should be " << countcurvednodes() << ". Have you correctly defined the node list?" << std::endl;
        abort();
    }
    if (quadrangleindex+1 > countquadrangularfaces())
    {
        std::cout << "Error: trying to get the " << quadrangleindex+1 << "th quadrangle in a " << gettypename() << " but there are only " << countquadrangularfaces() << std::endl;
        abort();
    }

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

std::vector<int> element::getstandardorientationreordering(void)
{    
    switch (getcurvedtypenumber())
    {
        // Nothing changes for points:
        case 0:
            return {0};
        // Lines stay untouched to avoid e.g. flip problems for the normal in 2D:
        case 1:
            return {0, 1};
        // For triangles and quadrangles make sure to circle around in the same 
        // direction to avoid e.g. flip problems for the normal:
        case 2:
        {
            // Get the index of the min node number in the straight element:
            int minnodenumindex = std::distance(curvednodelist.begin(), std::min_element(curvednodelist.begin(), curvednodelist.end()));
            return {minnodenumindex, (minnodenumindex+1)%3, (minnodenumindex+2)%3};
        }
        case 3:
        {
            // Get the index of the min node number in the straight element:
            int minnodenumindex = std::distance(curvednodelist.begin(), std::min_element(curvednodelist.begin(), curvednodelist.end()));
            return {minnodenumindex, (minnodenumindex+1)%4, (minnodenumindex+2)%4, (minnodenumindex+3)%4};
        }
        // For tetrahedra this is trivial:
        case 4:
        {
            // Get the indexes corresponding to the node numbers sorted ascendingly:
            std::vector<int> reorderingvector;
            myalgorithm::stablesort(curvednodelist, reorderingvector);
            return reorderingvector;
        }
        // For hexahedra make sure the numbering scheme is respected:
        case 5:
        {
            // Get the index of the min node number in the straight element:
            int minnodenumindex = std::distance(curvednodelist.begin(), std::min_element(curvednodelist.begin(), curvednodelist.end()));
    
            // Create a vector whose ith entry gives the indexes of the nodes connected to the ith node:
            std::vector<std::vector<int>> connectednodes = { {1,3,4},{2,0,5},{3,1,6},{0,2,7},{5,7,0},{6,4,1},{7,5,2},{4,6,3} };

            std::vector<int> output(8);
            // The first node is the smallest one:
            output[0] = minnodenumindex;
            // The second and third nodes are the smallest connected ones - but not the smallest node itself!
            for (int i = 1; i < 3; i++)
            {
                int minneighbourindex = connectednodes[output[i-1]][0];
                if (minneighbourindex == output[0])
                    minneighbourindex = connectednodes[output[i-1]][1];
                for (int j = 1; j < 3; j++)
                {
                    if (curvednodelist[connectednodes[output[i-1]][j]] < curvednodelist[minneighbourindex] && connectednodes[output[i-1]][j] != output[0])
                        minneighbourindex = connectednodes[output[i-1]][j];
                }
                output[i] = minneighbourindex;
            }
            // Set true for all nodes in the front face:
            std::vector<bool> isinfrontface(8,false);
            for (int i = 0; i < 3; i++)
                isinfrontface[output[i]] = true;
            // The fourth node is the last one in the current quad face.
            // For that it must be connected to the first and third nodes
            // and cannot be a previous node.
            for (int i = 0; i < 3; i++)
            {
                if (isinfrontface[connectednodes[output[0]][i]] == false)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        if (connectednodes[output[0]][i] == connectednodes[output[2]][j])
                            output[3] = connectednodes[output[0]][i];
                    }
                }
            }
            isinfrontface[output[3]] = true;
            // The last four nodes are on the opposite quad face:
            for (int i = 4; i < 8; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (isinfrontface[connectednodes[output[i-4]][j]] == false)
                        output[i] = connectednodes[output[i-4]][j];
                }
            }
            return output;
        }
        // For prisms:
        case 6:
        {
            // Get the index of the min node number in the straight element:
            int minnodenumindex = std::distance(curvednodelist.begin(), std::min_element(curvednodelist.begin(), curvednodelist.end()));
            // Give a 3 offset if the back face has the smallest node number:
            int offset = 0; 
            if (minnodenumindex > 2)
                offset = 3;
            
            std::vector<int> output(6);
            // Reorient the front or the back face:
            std::vector<int> reorderingvector;
            std::vector<int> trianglenodes = {curvednodelist[0+offset], curvednodelist[1+offset], curvednodelist[2+offset]};
            myalgorithm::stablesort(trianglenodes, reorderingvector);
            for (int i = 0; i < 3; i++)
                output[i] = reorderingvector[i]+offset; 
            for (int i = 3; i < 6; i++)
                output[i] = reorderingvector[i-3]+(offset+3)%6; 
            return output;
        }
        // For pyramids:
        case 7:
        {
            // Get the index of the min node number in the straight quad face:
            int minnodenumindex = std::distance(curvednodelist.begin(), std::min_element(curvednodelist.begin(), curvednodelist.begin()+4));
    
            if (curvednodelist[(minnodenumindex+1)%4] < curvednodelist[(minnodenumindex+3)%4])
                return {minnodenumindex, (minnodenumindex+1)%4, (minnodenumindex+2)%4, (minnodenumindex+3)%4, 4};
            else
                return {minnodenumindex, (minnodenumindex+3)%4, (minnodenumindex+2)%4, (minnodenumindex+1)%4, 4};
        }
        default:
            std::cout << "Error in 'element' object: element reordering is only defined for straight elements" << std::endl;
            abort;
    }
}

std::vector<int> element::getedgesreordering(std::vector<int> nodereordering)
{    
    std::vector<int> edges = getedgesdefinitionsbasedonnodes();
    std::vector<int> edgesreordered(edges.size()/2);
    
    // Loop on all reordered edges:
    for (int i = 0; i < edgesreordered.size(); i++)
    {
        std::vector<int> reorderededge = {nodereordering[edges[2*i]], nodereordering[edges[2*i+1]]};
        // Loop on all unreordered edges:
        for (int j = 0; j < edges.size()/2; j++)
        {
            std::vector<int> unreorderededge = {edges[2*j], edges[2*j+1]};
            if (reorderededge == unreorderededge || reorderededge[0] == unreorderededge[1] && reorderededge[1] == unreorderededge[0])
            {
                edgesreordered[i] = j;
                break;
            }
        }
    }
    return edgesreordered;
}

std::vector<int> element::gettriangularfacesreordering(std::vector<int> nodereordering)
{    
    int numtriangularfaces = counttriangularfaces();

    std::vector<int> faces = getfacesdefinitionsbasedonnodes();
    std::vector<int> facesreordered(numtriangularfaces);
    
    // Loop on all reordered faces:
    for (int i = 0; i < numtriangularfaces; i++)
    {
        std::vector<int> reorderedface = {nodereordering[faces[3*i]], nodereordering[faces[3*i+1]], nodereordering[faces[3*i+2]]};
        std::sort(reorderedface.begin(), reorderedface.end());
        // Loop on all unreordered faces:
        for (int j = 0; j < numtriangularfaces; j++)
        {
            std::vector<int> unreorderedface = {faces[3*j], faces[3*j+1], faces[3*j+2]};
            std::sort(unreorderedface.begin(), unreorderedface.end());
            if (reorderedface == unreorderedface)
            {
                facesreordered[i] = j;
                break;
            }
        }
    }
    return facesreordered;
}

std::vector<int> element::getquadrangularfacesreordering(std::vector<int> nodereordering)
{    
    int numtriangularfaces = counttriangularfaces();
    int numquadrangularfaces = countquadrangularfaces();

    std::vector<int> faces = getfacesdefinitionsbasedonnodes();
    // Remove any leading triangular face:
    for (int i = 0; i < numquadrangularfaces*4; i++)
        faces[i] = faces[i+numtriangularfaces*3];
    faces.resize(numquadrangularfaces*4);
    std::vector<int> facesreordered(numquadrangularfaces);
    
    // Loop on all reordered faces:
    for (int i = 0; i < numquadrangularfaces; i++)
    {
        std::vector<int> reorderedface = {nodereordering[faces[4*i]], nodereordering[faces[4*i+1]], nodereordering[faces[4*i+2]], nodereordering[faces[4*i+3]]};
        std::sort(reorderedface.begin(), reorderedface.end());
        // Loop on all unreordered faces:
        for (int j = 0; j < numquadrangularfaces; j++)
        {
            std::vector<int> unreorderedface = {faces[4*j], faces[4*j+1], faces[4*j+2], faces[4*j+3]};
            std::sort(unreorderedface.begin(), unreorderedface.end());
            if (reorderedface == unreorderedface)
            {
                facesreordered[i] = j;
                break;
            }
        }
    }
    return facesreordered;
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
                std::cout << "Error in 'element' object: coordinates of order 2 and above not defined for pyramids" << std::endl;
                abort();
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
}

std::vector<double> element::calculatecoordinates(std::vector<double>& refcoords, std::vector<double>& nodecoords)
{
    if (mypolynomials.count() > 0)
    {
        int numnodes = nodecoords.size()/3;
        int numrefs = refcoords.size()/3;
        
        std::vector<double> sf, evaluationpoint;
        std::vector<double> output(3*numrefs, 0.0);        

        for (int i = 0; i < numrefs; i++)
        {
            evaluationpoint = {refcoords[3*i+0],refcoords[3*i+1],refcoords[3*i+2]};
            
            mypolynomials.evalatsingle(evaluationpoint, sf);

            for (int c = 0; c < numnodes; c++)
            {
                output[3*i+0] += nodecoords[3*c+0] * sf[c];
                output[3*i+1] += nodecoords[3*c+1] * sf[c];
                output[3*i+2] += nodecoords[3*c+2] * sf[c];
            }
        }
        return output;
    }
    else
    {
        lagrangeformfunction lff(gettypenumber(), getcurvatureorder(), {});
        mypolynomials = polynomials(lff.getformfunctionpolynomials());
        return calculatecoordinates(refcoords, nodecoords);
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
        std::vector<double> oncesplitcoords = cursplitcoords[gettypenumber()];
        fullsplit(n-1, splitcoords, oncesplitcoords);
        // Treat the tetrahedra from the split pyramids:
        if (gettypenumber() == 7)
        {
            std::vector<std::vector<double>> tetsplitcoords;
            fullsplit(n-1, tetsplitcoords, cursplitcoords[4]);
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
}

void element::write(std::string filename, std::vector<double> coords)
{
    int ncn = countcurvednodes();
    int numelems = coords.size()/3/ncn;
    
    densematrix xcoords(numelems, ncn);
    densematrix ycoords(numelems, ncn);
    densematrix zcoords(numelems, ncn);
    
    densematrix vals(numelems, ncn);

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

