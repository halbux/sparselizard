#include "gmshinterface.h"


void gmshinterface::readfromfile(std::string name, nodes& mynodes, elements& myelements, physicalregions& myphysicalregions)
{	
	std::string currentline;

	// 'file' cannot take a std::string argument --> name.c_str():
	std::ifstream meshfile (name.c_str());
	if (meshfile.is_open())
	{
		// Move to the mesh version section:
		double formatversion;
		while (std::getline(meshfile, currentline))
		{
			if (currentline == "$MeshFormat")
			{
				std::getline(meshfile, currentline);
				// Get version number:
				mystring stringobject(currentline);
				formatversion = std::stod(stringobject.getstringtonextwhitespace());
				break;
			}
		}
		// Give an error if version is not supported:
		if (formatversion >= 4)
		{
			std::cout << "Error in 'gmshinterface': GMSH format " << formatversion << " is not supported yet (might still be experimental)." << std::endl;
			std::cout << "Use GMSH 3 instead or the built-in mesher." << std::endl;
			abort();
		}

		// Move to the node section and read the number of nodes:
		int numberofnodes;
		while (std::getline(meshfile, currentline))
		{
			if (currentline == "$Nodes")
			{
				std::getline(meshfile, currentline);
				numberofnodes = std::stoi(currentline);
				break;
			}
		}

		// We are in the node section and now get the node coordinates (doubles):
        mynodes.setnumber(numberofnodes);
		std::vector<double>* nodecoordinates = mynodes.getcoordinates();
		
		for (int i = 0; i < numberofnodes; i++)
		{
			std::getline(meshfile, currentline);
            mystring stringobject(currentline);
			
			// The first number in the line is an integer (the node number). We skip it:
			stringobject.jumptonextwhitespace();
			// Read now the x, y and z node coordinate:
			nodecoordinates->at(3*i+0) = std::stod(stringobject.getstringtonextwhitespace());
			nodecoordinates->at(3*i+1) = std::stod(stringobject.getstringtonextwhitespace());
			nodecoordinates->at(3*i+2) = std::stod(stringobject.getstringtonextwhitespace());
		}

		// Move to the element section and read the number of elements:
		int numberofelements;
		while (std::getline(meshfile, currentline))
		{
			if (currentline == "$Elements")
			{
				std::getline(meshfile, currentline);
				numberofelements = std::stoi(currentline);
				break;
			}
		}

		// We are in the element section and now fill in 'myelements' and 'myphysicalregions':
        int elementindexincurrenttype = 0;
		std::vector<int> nodesinpreviouselement = {};

        for (int i = 0; i < numberofelements; i++)
		{
			std::getline(meshfile, currentline);
			mystring stringobject(currentline);
			
			// The first number in the line is an integer (the element number). We skip it:
			stringobject.jumptonextwhitespace();
			// Read now the element type number and define the element object:
			int currentcurvedelementtype = gmshinterface::convertgmshelementtypenumber(std::stoi(stringobject.getstringtonextwhitespace()));
            element elementobject(currentcurvedelementtype);
            // Get the uncurved element type number:
            int currentelementtype = elementobject.gettypenumber();
            int curvatureorder = elementobject.getcurvatureorder();
            
			// Read now the number of parameters (integers). First parameter is the physical region, we skip the other ones:
			int numberofparameters = std::stoi(stringobject.getstringtonextwhitespace());
			// Read only the first parameter:
            int currentphysicalregionnumber;
			for (int j = 0; j < numberofparameters; j++)
			{
				// Get the current physical region number:
				if (j == 0)
					currentphysicalregionnumber = std::stoi(stringobject.getstringtonextwhitespace());
				else
					stringobject.jumptonextwhitespace();
			}
            // Get the physical region object associated to 'currentphysicalregionnumber':
            physicalregion* currentphysicalregion = myphysicalregions.get(currentphysicalregionnumber);
			// Read now the node number list in the element. The number of nodes depends on the element:
			std::vector<int>  nodesincurrentelement(elementobject.countcurvednodes());
			for (int j = 0; j < elementobject.countcurvednodes(); j++)
				nodesincurrentelement[j] = std::stoi(stringobject.getstringtonextwhitespace()) - 1; // -1 to start numbering nodes at 0
            
			// If the current element is the same as the previous one then the element is already
			// in the 'myelements' object and it only has to be added to the physical region. 
			// This can often occurs in some .msh files.
			if (nodesinpreviouselement == nodesincurrentelement)
			{
                currentphysicalregion->addelement(currentelementtype, elementindexincurrenttype);
				continue;
			}
			
			// In case it is a new element add the element AND its physical region.
            elementindexincurrenttype = myelements.add(currentelementtype, curvatureorder, nodesincurrentelement);
            currentphysicalregion->addelement(currentelementtype, elementindexincurrenttype);

			nodesinpreviouselement = nodesincurrentelement;
		}
		
		meshfile.close();
	}
	else 
	{
		std::cout << "Unable to open file " << name << " or file not found" << std::endl;
		abort();
	}
}

void gmshinterface::writetofile(std::string name, nodes& mynodes, elements& myelements, physicalregions& myphysicalregions, disjointregions& mydisjointregions)
{	
	// 'file' cannot take a std::string argument --> name.c_str():
	std::ofstream meshfile (name.c_str());
	if (meshfile.is_open())
	{
		// To write all doubles with enough digits to the file:
		meshfile << std::setprecision(17);
		
		// Write the header:
		meshfile << "$MeshFormat\n";
		meshfile << "2.2 0 8\n";
		meshfile << "$EndMeshFormat\n";
		
		// Write the node section:
		meshfile << "$Nodes\n";
		meshfile << mynodes.count() << "\n";
		// Write the node coordinates:		
        std::vector<double>* nodecoordinates = mynodes.getcoordinates();
        
		for (int i = 0; i < mynodes.count(); i++)
			meshfile << i+1 << " " << nodecoordinates->at(3*i+0) << " " << nodecoordinates->at(3*i+1) << " " << nodecoordinates->at(3*i+2) << "\n";
		meshfile << "$EndNodes\n";
		
		// Write the element section:
		meshfile << "$Elements\n";
		
		// Write the number of elements.
		// The number of elements is equal to the total number of elements in the physical regions.
		// Elements that are in several physical regions at the same time are counted 
		// multiple times, which is ok in the .msh format.
        int numberofelements = myphysicalregions.countelements();
		meshfile << numberofelements << "\n";
		
		// Write all element lines to the file.
		int elementnumberinfile = 1;
		// Iterate through all physical regions:
		for (int i = 0; i < myphysicalregions.count(); i++)
		{
            // Get the physical region number corresponding to index i:
            int physicalregionnumber = myphysicalregions.getnumber(i);
            // Get the physical region object:
            physicalregion* physicalregionobject = myphysicalregions.getatindex(i);
            // Get all disjoint regions inside the physical region:
            std::vector<int> alldisjointregions = physicalregionobject->getdisjointregions();

			// Iterate on all disjoint regions in the physical region:
			for (int h = 0; h < alldisjointregions.size(); h++)
			{
                // Get the unique element type in the disjoint region:
                int typenumber = mydisjointregions.getelementtypenumber(alldisjointregions[h]);
                // Get the corresponding high order element type number:
                element myelement(typenumber, myelements.getcurvatureorder());
                int curvedtypenumber = myelement.getcurvedtypenumber();
                int numberofcurvednodes = myelement.countcurvednodes();
                
                // Get the range begin and end in the disjoint region:
                int rangebegin = mydisjointregions.getrangebegin(alldisjointregions[h]);
                int rangeend = mydisjointregions.getrangeend(alldisjointregions[h]);

                // Iterate on all elements:
				for (int i = rangebegin; i <= rangeend; i++)
				{
                    // Write [element number in file, type number, number of parameters (2 here), physical region number, physical region number]:
					// Note: physical region number appears twice. If a single parameter is provided the
					// physical regions will not be displayable separately when the .msh file is opened in GMSH.
					meshfile << elementnumberinfile << " " << gmshinterface::converttogmshelementtypenumber(curvedtypenumber) << " 2 " << physicalregionnumber << " " << physicalregionnumber;
					// Write all nodes in the element.
                    for (int nodeindex = 0; nodeindex < numberofcurvednodes; nodeindex++)
                        meshfile << " " << myelements.getsubelement(0, typenumber, i, nodeindex) + 1; // +1 to start numbering nodes at 1
                    
					elementnumberinfile = elementnumberinfile + 1;
					meshfile << "\n";
				}
			}
		}
			
		meshfile << "$EndElements";
		
		meshfile.close();
	}
	else 
	{
		std::cout << "Unable to write to file " << name << " or file not found" << std::endl;
		abort();
	}
}

void gmshinterface::openview(std::string name, std::string viewname, double timetag, bool overwrite)
{	
	// 'file' cannot take a std::string argument --> name.c_str():
    std::ofstream posfile;
    if (overwrite)
        posfile.open(name.c_str());
    else
        posfile.open(name.c_str(), std::ios::out | std::ios::app );
        
	if (posfile.is_open())
	{
		// To write all doubles with enough digits to the file:
		posfile << std::setprecision(17);

		// Write the header:
        posfile << "View \"" << viewname << "\" {\n";
        // Write the time tag:
        posfile << "TIME{" << timetag << "};\n";
        
		posfile.close();
	}
	else 
	{
		std::cout << "Unable to write to file " << name << " or file not found" << std::endl;
		abort();
	}
}

void gmshinterface::appendtoview(std::string name, int elementtypenumber, densematrix coordx, densematrix coordy, densematrix coordz, densematrix compxinterpolated)
{	
	// 'file' cannot take a std::string argument --> name.c_str():
	std::ofstream posfile;
    posfile.open(name.c_str(), std::ios::out | std::ios::app );
	if (posfile.is_open())
	{
		// To write all doubles with enough digits to the file:
		posfile << std::setprecision(17);
		
        // Get the character identifying the element:
        char elementidentifier = getelementidentifierinposformat(elementtypenumber);
        
        // Write all elements:
        for (int elem = 0; elem < compxinterpolated.countrows(); elem++)
        {
            posfile << 'S' << elementidentifier << "(";
            // Write the coordinates of every evaluation point:
            for (int i = 0; i < coordx.countcolumns(); i++)
            {
                posfile << coordx.getvalue(elem,i) << "," << coordy.getvalue(elem,i) << "," << coordz.getvalue(elem,i);
                if (i < coordx.countcolumns()-1)
                    posfile << ",";
            }
            // Write the field value at every evaluation point:
            posfile << ")\n{";
            for (int i = 0; i < compxinterpolated.countcolumns(); i++)
            {
                posfile << compxinterpolated.getvalue(elem,i);
                if (i < compxinterpolated.countcolumns()-1)
                    posfile << ",";
            }
            posfile << "};\n";
        }
        
		posfile.close();
	}
	else 
	{
		std::cout << "Unable to write to file " << name << " or file not found" << std::endl;
		abort();
	}
}

void gmshinterface::appendtoview(std::string name, int elementtypenumber, densematrix coordx, densematrix coordy, densematrix coordz, densematrix compxinterpolated, densematrix compyinterpolated, densematrix compzinterpolated)
{	
	// 'file' cannot take a std::string argument --> name.c_str():
	std::ofstream posfile;
    posfile.open(name.c_str(), std::ios::out | std::ios::app );
	if (posfile.is_open())
	{
		// To write all doubles with enough digits to the file:
		posfile << std::setprecision(17);
		
        // Get the character identifying the element:
        char elementidentifier = getelementidentifierinposformat(elementtypenumber);
        
        // Write all elements:
        for (int elem = 0; elem < compxinterpolated.countrows(); elem++)
        {
            posfile << 'V' << elementidentifier << "(";
            // Write the coordinates of every evaluation point:
            for (int i = 0; i < coordx.countcolumns(); i++)
            {
                posfile << coordx.getvalue(elem,i) << "," << coordy.getvalue(elem,i) << "," << coordz.getvalue(elem,i);
                if (i < coordx.countcolumns()-1)
                    posfile << ",";
                
            }
            // Write the field value at every evaluation point:
            posfile << ")\n{";
            for (int i = 0; i < compxinterpolated.countcolumns(); i++)
            {
                posfile << compxinterpolated.getvalue(elem,i) << "," << compyinterpolated.getvalue(elem,i) << "," << compzinterpolated.getvalue(elem,i);
                if (i < compxinterpolated.countcolumns()-1)
                    posfile << ",";
            }
            posfile << "};\n";
        }
        
		posfile.close();
	}
	else 
	{
		std::cout << "Unable to write to file " << name << " or file not found" << std::endl;
		abort();
	}
}

void gmshinterface::writeinterpolationscheme(std::string name, std::vector<std::vector<polynomial>> poly)
{	
	// 'file' cannot take a std::string argument --> name.c_str():
	std::ofstream posfile;
    posfile.open(name.c_str(), std::ios::out | std::ios::app );
	if (posfile.is_open())
	{
        posfile << "\nINTERPOLATION_SCHEME";
        
        for (int m = 0; m < poly.size(); m++)
        {
            posfile << "\n{\n";
            
            // Print the polynomial coefficients:
            for (int p = 0; p < poly[m].size(); p++)
            {
                posfile << "  {";

                std::vector<std::vector<std::vector<double>>> polyformfunctions = poly[m][p].get();

                for (int i = 0; i < polyformfunctions.size(); i++)
                {
                    for (int j = 0; j < polyformfunctions[i].size(); j++)
                    {
                        for (int k = 0; k < polyformfunctions[i][j].size(); k++)
                        {
                            if (i == polyformfunctions.size() - 1 && j == polyformfunctions[i].size() - 1 && k == polyformfunctions[i][j].size() - 1)
                                posfile << polyformfunctions[i][j][k];
                            else
                                posfile << polyformfunctions[i][j][k] << ",";
                        }
                    }
                }

                if (p == poly[m].size() - 1)
                    posfile << "}\n";
                else
                    posfile << "},\n";
            }

            posfile << "}\n{\n";

            // Print the list of monomials:
            std::vector<std::vector<std::vector<double>>> polyformfunctions = poly[m][0].get();

            for (int i = 0; i < polyformfunctions.size(); i++)
            {
                for (int j = 0; j < polyformfunctions[i].size(); j++)
                {
                    for (int k = 0; k < polyformfunctions[i][j].size(); k++)
                    {
                        posfile << "  {" << i << "," << j << "," << k << "}";
                        if (i == polyformfunctions.size() - 1 && j == polyformfunctions[i].size() - 1 && k == polyformfunctions[i][j].size() - 1)
                            posfile << "\n";
                        else
                            posfile << ",\n";
                    }
                }
            }
            posfile << "}";
        }

        posfile << ";\n\n";

		posfile.close();
	}
	else 
	{
		std::cout << "Unable to write to file " << name << " or file not found" << std::endl;
		abort();
	}
}

void gmshinterface::closeview(std::string name)
{	
	// 'file' cannot take a std::string argument --> name.c_str():
	std::ofstream posfile;
    posfile.open(name.c_str(), std::ios::out | std::ios::app );
	if (posfile.is_open())
	{
        posfile << "};\n\n";
        
		posfile.close();
	}
	else 
	{
		std::cout << "Unable to write to file " << name << " or file not found" << std::endl;
		abort();
	}
}

int gmshinterface::convertgmshelementtypenumber(int gmshtypenumber)
{	
	// Elements 1 to 14 are identical:
	if (gmshtypenumber > 0 && gmshtypenumber < 15)
		return gmshtypenumber;
	else
	{
		switch (gmshtypenumber)
		{
			// Point:
			case 15:
				return 0;
                
			// Line order 3:
			case 26:
				return 15;
			// Line order 4:
			case 27:
				return 22;
			// Line order 5:
			case 28:
				return 29;
                
			// Triangle order 3:
			case 21:
				return 16;
			// Triangle order 4:
			case 23:
				return 23;
			// Triangle order 5:
			case 25:
				return 30;
                
			// Quadrangle order 3:
			case 36:
				return 17;
			// Quadrangle order 4:
			case 37:
				return 24;
			// Quadrangle order 5:
			case 38:
				return 31;
                
			// Tetrahedron order 3:
			case 29:
				return 18;
			// Tetrahedron order 4:
			case 30:
				return 25;
			// Tetrahedron order 5:
			case 31:
				return 32;
                
			// Hexahedron order 3:
			case 92:
				return 19;
			// Hexahedron order 4:
			case 93:
				return 26;
			// Hexahedron order 5:
			case 94:
				return 33;
                
			// Prism order 3:
			case 90:
				return 20;
                
			default:
				std::cout << "Error in 'gmshinterface' namespace: trying to use a GMSH element that is undefined in this code." << std::endl;
				abort();
		}
	}
}

int gmshinterface::converttogmshelementtypenumber(int ourtypenumber)
{
	// Elements 1 to 14 are identical:
	if (ourtypenumber > 0 && ourtypenumber < 15)
		return ourtypenumber;
	else
	{
		switch (ourtypenumber)
		{
			// Point:
			case 0:
				return 15;
                
			// Line order 3:
			case 15:
				return 26;
			// Line order 4:
			case 22:
				return 27;
			// Line order 5:
			case 29:
				return 28;
                
			// Triangle order 3:
			case 16:
				return 21;
			// Triangle order 4:
			case 23:
				return 23;
			// Triangle order 5:
			case 30:
				return 25;
                
			// Quadrangle order 3:
			case 17:
				return 36;
			// Quadrangle order 4:
			case 24:
				return 37;
			// Quadrangle order 5:
			case 31:
				return 38;
                
			// Tetrahedron order 3:
			case 18:
				return 29;
			// Tetrahedron order 4:
			case 25:
				return 30;
			// Tetrahedron order 5:
			case 32:
				return 31;
                
			// Hexahedron order 3:
			case 19:
				return 92;
			// Hexahedron order 4:
			case 26:
				return 93;
			// Hexahedron order 5:
			case 33:
				return 94;
                
			// Prism order 3:
			case 20:
				return 90;
                
			default:
				std::cout << "Error in 'gmshinterface' namespace: trying to use a GMSH element that is undefined in this code." << std::endl;
				abort();
		}
    }
}

char gmshinterface::getelementidentifierinposformat(int ourtypenumber)
{	
    element myelement(ourtypenumber);
    
    // Switch on the straight type number:
	switch (myelement.gettypenumber())
	{
		// Point:
		case 0:
			return 'P';
		// Line:
		case 1:
			return 'L';
		// Triangle:
		case 2:
			return 'T';
		// Quadrangle:
		case 3:
			return 'Q';
		// Tetrahedron:
		case 4:
			return 'S';
		// Hexahedron:
		case 5:
			return 'H';
		// Prism:
		case 6:
			return 'I';
		// Pyramid:
		case 7:
			return 'Y';
	}
}

