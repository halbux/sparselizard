#include "orientation.h"

int orientation::gettotalorientation(int elementtypenumber, std::vector<int>& nodelist)
{
    int totalorientation = 0;
    int factor = 1;
    
    // Get the orientation of all edges and faces in the element:
    std::vector<int> orientationofedges = getedgesorientationsinelement(elementtypenumber, nodelist);
    std::vector<int> orientationoffaces = getfacesorientationsinelement(elementtypenumber, nodelist);
    
    // Create the total orientation number:
    for (int i = 0; i < orientationofedges.size(); i++)
    {
        totalorientation = totalorientation + orientationofedges[i]*factor;
        factor = factor * 2;
    }
    for (int i = 0; i < orientationoffaces.size(); i++)
    {
        totalorientation = totalorientation + orientationoffaces[i]*factor;
        factor = factor * 8;
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

std::vector<int> orientation::getedgesorientationsinelement(int elementtypenumber, std::vector<int>& curvednodelist)
{
    element myelement(elementtypenumber);
    
    int numberofedges = myelement.countedges();
    std::vector<int> alledgesorientations(numberofedges);
    std::vector<int> edgesdefinitionbasedonnodes = myelement.getedgesdefinitionsbasedonnodes();
    std::vector<int> nodesinedges(2);
    
    // Loop on all edges:
    for (int i = 0; i < numberofedges; i++)
    {
        // Create the edge node list:
        nodesinedges[0] = curvednodelist[edgesdefinitionbasedonnodes[2*i+0]];
        nodesinedges[1] = curvednodelist[edgesdefinitionbasedonnodes[2*i+1]];
        
        alledgesorientations[i] = getorientationofedge(nodesinedges);
    }
    return alledgesorientations;
}

std::vector<int> orientation::getfacesorientationsinelement(int elementtypenumber, std::vector<int>& curvednodelist)
{
    element myelement(elementtypenumber);
    
    int numberoffaces = myelement.countfaces();
    int numberoftriangularfaces = myelement.counttriangularfaces();
    
    std::vector<int> allfacesorientations(numberoffaces);
    std::vector<int> facesdefinitionbasedonnodes = myelement.getfacesdefinitionsbasedonnodes();
    
    // Loop on all faces:
    for (int i = 0; i < numberoffaces; i++)
    {
        if (myelement.istriangularface(i))
        {
            std::vector<int> nodesinface(3);

            // Create the triangle node list:
            nodesinface[0] = curvednodelist[facesdefinitionbasedonnodes[3*i+0]];
            nodesinface[1] = curvednodelist[facesdefinitionbasedonnodes[3*i+1]];
            nodesinface[2] = curvednodelist[facesdefinitionbasedonnodes[3*i+2]];
        
            allfacesorientations[i] = getorientationoftriangle(nodesinface);
        }
        else
        {
            std::vector<int> nodesinface(4);
            int offset = 3*numberoftriangularfaces;

            // Create the quadrangle node list:
            nodesinface[0] = curvednodelist[facesdefinitionbasedonnodes[offset+4*(i-numberoftriangularfaces)+0]];
            nodesinface[1] = curvednodelist[facesdefinitionbasedonnodes[offset+4*(i-numberoftriangularfaces)+1]];
            nodesinface[2] = curvednodelist[facesdefinitionbasedonnodes[offset+4*(i-numberoftriangularfaces)+2]];
            nodesinface[3] = curvednodelist[facesdefinitionbasedonnodes[offset+4*(i-numberoftriangularfaces)+3]];

            allfacesorientations[i] = getorientationofquadrangle(nodesinface);
        }
    }
    return allfacesorientations;
}

int orientation::getorientationofedge(std::vector<int>& physicalnodesinedge)
{
    if (physicalnodesinedge[0] > physicalnodesinedge[1])
        return 0;
    else
        return 1;
}

int orientation::getorientationoftriangle(std::vector<int>& physicalnodesintriangle)
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

int orientation::getorientationofquadrangle(std::vector<int>& physicalnodesinquadrangle)
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

