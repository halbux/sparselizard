#include "element.h"


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
		std::cout << "Element type number is " << number << ", curvature order is " << curvatureorder << std::endl;
		abort();
	}
	if (curvatureorder <= 0)
	{
		std::cout << "Error in 'element' object: can not define a negative or 0 curvature order" << std::endl;
		std::cout << "Element type number is " << number << ", curvature order is " << curvatureorder << std::endl;
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
		return curvedtypenumber/7.0;
	else
		return (curvedtypenumber - curvedtypenumber%7)/7.0 + 1;		
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
			curvednumberofnodes = 0.5*( pow(order+1,2) + (order+1) );
			break;
		// Quadrangle:
		case 3:
			curvednumberofnodes = pow(order + 1,2);
			break;
		// Tetrahedron:
		case 4:
			curvednumberofnodes = 1.0/6.0*(order+1)*(order+2)*(order+3);
            break;
		// Hexahedron:
		case 5:
			curvednumberofnodes = pow(order + 1,3);
			break;
		// Prism:
		case 6:
			curvednumberofnodes = 0.5*(order + 1)*( ( pow(order+1,2) + (order+1) ) );
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

int element::countnodes(void)
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

int element::countedges(void)
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

int element::countfaces(void)
{
	return counttriangularfaces() + countquadrangularfaces();
}

int element::counttriangularfaces(void)
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

int element::countquadrangularfaces(void)
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

int element::countvolumes(void)
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
			return 0;
		// Tetrahedron:
		case 4:
			return 1;
		// Hexahedron:
		case 5:
			return 1;
		// Prism:
		case 6:
			return 1;
		// Pyramid:
		case 7:
			return 1;
	}
}

int element::count(int dim)
{
    switch (dim)
    {
		case 0:
			return countnodes();
		case 1:
			return countedges();
		case 2:
			return countfaces();
		case 3:
			return countvolumes();
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

