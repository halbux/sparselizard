#include "orientation.h"

std::vector<int> orientation::gettotalorientation(int elementtypenumber, std::vector<long long int>& nodelist)
{
    element myelement(elementtypenumber);
    int numberofnodes = myelement.countnodes();
    int numberofedges = myelement.countedges();
    int numberoffaces = myelement.countfaces();
    int numberofelements = nodelist.size()/numberofnodes;
    
    std::vector<int> totalorientation(numberofelements, 0);
    
    // Get the orientation of all edges and faces in the element:
    std::vector<int> orientationofedges = getedgesorientationsinelement(elementtypenumber, nodelist);
    std::vector<int> orientationoffaces = getfacesorientationsinelement(elementtypenumber, nodelist);
    
    // Create the total orientation number:
    for (int e = 0; e < numberofelements; e++)
    {
        int factor = 1;
        for (int i = 0; i < numberofedges; i++)
        {
            totalorientation[e] += orientationofedges[e*numberofedges+i]*factor;
            factor *= 2;
        }
        for (int i = 0; i < numberoffaces; i++)
        {
            totalorientation[e] += orientationoffaces[e*numberoffaces+i]*factor;
            factor *= 8;
        }
    }
    
    return totalorientation;
}

int orientation::countorientations(int elemtypenum)
{
    switch (elemtypenum)
    {
        case 0:
            return 1;
        case 1:
            return 2;
        case 2:
            return 6;
        case 3:
            return 8;
        case 4:
            return 1;
        case 5:
            return 1;
        case 6:
            return 1;
        case 7:
            return 1;
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<int> orientation::getedgesorientationsfromtotalorientation(int totalorientation, int elementtypenumber)
{
    element myelement(elementtypenumber);
    int numberofedges = myelement.countedges();

    std::vector<int> edgesorientations(numberofedges);

    for (int edge = 0; edge < numberofedges; edge++)
    {
        int remainder = totalorientation%2;
        edgesorientations[edge] = remainder;    
        totalorientation = (totalorientation-remainder)/2;
    }
        
    return edgesorientations;
}

std::vector<int> orientation::getfacesorientationsfromtotalorientation(int totalorientation, int elementtypenumber)
{
    element myelement(elementtypenumber);
    int numberofedges = myelement.countedges();
    int numberoffaces = myelement.countfaces();
    
    int twotothenumberofedges = pow(2,numberofedges);
    
    // Remove the part that corresponds to the edges orientations:
    totalorientation = (totalorientation-totalorientation%(twotothenumberofedges)) / twotothenumberofedges;

    std::vector<int> facesorientations(numberoffaces);

    for (int face = 0; face < numberoffaces; face++)
    {
        int remainder = totalorientation%8;
        facesorientations[face] = remainder;    
        totalorientation = (totalorientation-remainder)/8;
    }
        
    return facesorientations;
}

std::vector<std::vector<int>> orientation::getreorderingtoreferenceedgeorientation(void)
{
    std::vector<std::vector<int>> reordering = {{0,1},{1,0}};
    return reordering;
}

std::vector<std::vector<int>> orientation::getreorderingtoreferencetriangularfaceorientation(void)
{
    std::vector<std::vector<int>> reordering = {{0,1,2},{0,2,1},{1,2,0},{1,0,2},{2,0,1},{2,1,0}};
    return reordering;
}

std::vector<std::vector<int>> orientation::getreorderingtoreferencequadrangularfaceorientation(void)
{
    std::vector<std::vector<int>> reordering = {{0,1,2,3},{0,3,2,1},{1,2,3,0},{1,0,3,2},{2,3,0,1},{2,1,0,3},{3,0,1,2},{3,2,1,0}};
    return reordering;
}

std::vector<int> orientation::getedgesorientationsinelement(int elementtypenumber, std::vector<long long int>& nodelist)
{
    element myelement(elementtypenumber);
    int numberofnodes = myelement.countnodes();
    int numberofedges = myelement.countedges();
    int numberofelements = nodelist.size()/numberofnodes;
    
    std::vector<int> alledgesorientations(numberofelements*numberofedges);
    std::vector<int> edgesdefinitionbasedonnodes = myelement.getedgesdefinitionsbasedonnodes();
    std::vector<long long int> nodesinedges(2);
    
    // Loop on all elements:
    for (int e = 0; e < numberofelements; e++)
    {
        // Loop on all edges:
        for (int i = 0; i < numberofedges; i++)
        {
            // Create the edge node list:
            nodesinedges[0] = nodelist[e*numberofnodes+edgesdefinitionbasedonnodes[2*i+0]];
            nodesinedges[1] = nodelist[e*numberofnodes+edgesdefinitionbasedonnodes[2*i+1]];
            
            alledgesorientations[e*numberofedges+i] = getorientationofedge(nodesinedges);
        }
    }
    return alledgesorientations;
}

std::vector<int> orientation::getfacesorientationsinelement(int elementtypenumber, std::vector<long long int>& nodelist)
{
    element myelement(elementtypenumber);
    int numberofnodes = myelement.countnodes();
    int numberoffaces = myelement.countfaces();
    int numberoftriangularfaces = myelement.counttriangularfaces();
    int numberofelements = nodelist.size()/numberofnodes;
    
    std::vector<int> allfacesorientations(numberofelements*numberoffaces);
    std::vector<int> facesdefinitionbasedonnodes = myelement.getfacesdefinitionsbasedonnodes();
    
    // Loop on all elements:
    for (int e = 0; e < numberofelements; e++)
    {
        // Loop on all faces:
        for (int i = 0; i < numberoffaces; i++)
        {
            if (i < numberoftriangularfaces)
            {
                std::vector<long long int> nodesinface(3);

                // Create the triangle node list:
                nodesinface[0] = nodelist[e*numberofnodes+facesdefinitionbasedonnodes[3*i+0]];
                nodesinface[1] = nodelist[e*numberofnodes+facesdefinitionbasedonnodes[3*i+1]];
                nodesinface[2] = nodelist[e*numberofnodes+facesdefinitionbasedonnodes[3*i+2]];
            
                allfacesorientations[e*numberoffaces+i] = getorientationoftriangle(nodesinface);
            }
            else
            {
                std::vector<long long int> nodesinface(4);
                int offset = 3*numberoftriangularfaces;

                // Create the quadrangle node list:
                nodesinface[0] = nodelist[e*numberofnodes+facesdefinitionbasedonnodes[offset+4*(i-numberoftriangularfaces)+0]];
                nodesinface[1] = nodelist[e*numberofnodes+facesdefinitionbasedonnodes[offset+4*(i-numberoftriangularfaces)+1]];
                nodesinface[2] = nodelist[e*numberofnodes+facesdefinitionbasedonnodes[offset+4*(i-numberoftriangularfaces)+2]];
                nodesinface[3] = nodelist[e*numberofnodes+facesdefinitionbasedonnodes[offset+4*(i-numberoftriangularfaces)+3]];

                allfacesorientations[e*numberoffaces+i] = getorientationofquadrangle(nodesinface);
            }
        }
    }
    return allfacesorientations;
}

int orientation::getorientationofedge(std::vector<long long int>& physicalnodesinedge)
{
    if (physicalnodesinedge[0] > physicalnodesinedge[1])
        return 0;
    else
        return 1;
}

int orientation::getorientationoftriangle(std::vector<long long int>& physicalnodesintriangle)
{
    // Get the position of the maximum:
    int maxpos = 0;
    for (int i = 1; i < 3; i++)
    {
        if (physicalnodesintriangle[i] > physicalnodesintriangle[maxpos])
            maxpos = i;
    }
    if (physicalnodesintriangle[(maxpos+1)%3] > physicalnodesintriangle[(maxpos+2)%3])
        return 2*maxpos;
    else
        return 2*maxpos + 1;
}

int orientation::getorientationofquadrangle(std::vector<long long int>& physicalnodesinquadrangle)
{
    // Get the position of the maximum:
    int maxpos = 0;
    for (int i = 1; i < 4; i++)
    {
        if (physicalnodesinquadrangle[i] > physicalnodesinquadrangle[maxpos])
            maxpos = i;
    }
    if (physicalnodesinquadrangle[(maxpos+1)%4] > physicalnodesinquadrangle[(maxpos+3)%4])
        return 2*maxpos;
    else
        return 2*maxpos + 1;
}

