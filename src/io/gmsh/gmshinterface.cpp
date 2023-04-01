#include "gmshinterface.h"
#include "universe.h"


#ifndef HAVE_GMSH
void gmshinterface::readfromapi(nodes& mynodes, elements& myelements, physicalregions& myphysicalregions)
{    
    logs log;
    log.msg() << "Error in 'gmshinterface' namespace: GMSH API is not available" << std::endl;
    log.error();
}

void gmshinterface::readwithapi(std::string name, nodes& mynodes, elements& myelements, physicalregions& myphysicalregions)
{
    logs log;
    log.msg() << "Error in 'gmshinterface' namespace: GMSH API is not available" << std::endl;
    log.error();
}
#endif
#ifdef HAVE_GMSH
#include "gmsh.h"
void gmshinterface::readfromapi(nodes& mynodes, elements& myelements, physicalregions& myphysicalregions)
{    

    ///// Get all nodes and coordinates:
    
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords, parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoord, -1, -1, false, false);

    if (nodeTags.size() == 0)
    {
        logs log;
        log.msg() << "Error in 'gmshinterface' namespace: no mesh node found in mesh loaded from gmsh api" << std::endl;
        log.error();
    }

    int numberofnodes = 0;
    
    int maxnodetag = *std::max_element(nodeTags.begin(), nodeTags.end());
    std::vector<int> noderenumbering(maxnodetag+1, -1);
    
    // The coordinates follow the nodeTags order:
    std::vector<int> nodeTagsrenumbering(maxnodetag+1, -1);
    for (int i = 0; i < nodeTags.size(); i++)
        nodeTagsrenumbering[nodeTags[i]] = i;
    
        
    ///// Get all physical region numbers:
    
    gmsh::vectorpair dimTags;
    gmsh::model::getPhysicalGroups(dimTags, -1);

    int numphysregs = dimTags.size();
    if (numphysregs == 0)
    {
        logs log;
        log.msg() << "Error in 'gmshinterface' namespace: no physical region found in mesh loaded from gmsh api" << std::endl;
        log.error();
    }
    std::vector<int> allphysregsdims(numphysregs);
    std::vector<int> allphysregs(numphysregs);
    for (int i = 0; i < numphysregs; i++)
    {
        allphysregsdims[i] = dimTags[i].first;
        allphysregs[i] = dimTags[i].second;
    }
        
        
    ///// Get the entities in each physical region:
    
    std::vector<std::vector<int>> entitiesinphysreg(numphysregs);
    for (int i = 0; i < numphysregs; i++)
        gmsh::model::getEntitiesForPhysicalGroup(allphysregsdims[i], allphysregs[i], entitiesinphysreg[i]);
        
        
    ///// Get the elements in each physical region:
    
    for (int i = 0; i < numphysregs; i++)
    {
        int dim = allphysregsdims[i];
        int currentphysicalregionnumber = allphysregs[i];
        physicalregion* currentphysicalregion = myphysicalregions.get(universe::physregshift*(dim+1) + currentphysicalregionnumber);
        
        // Loop on all entities in the physical region:
        for (int j = 0; j < entitiesinphysreg[i].size(); j++)
        {
            std::vector<std::vector<std::size_t>> elemnodeTags, elementTags;
            std::vector<int> elementTypes;

            gmsh::model::mesh::getElements(elementTypes, elementTags, elemnodeTags, dim, entitiesinphysreg[i][j]);

            for (int t = 0; t < elementTypes.size(); t++)
            {
                int currentcurvedelementtype = convertgmshelementtypenumber(elementTypes[t]);
                element elementobject(currentcurvedelementtype);
                // Get the uncurved element type number:
                int currentelementtype = elementobject.gettypenumber();
                int curvatureorder = elementobject.getcurvatureorder();
                int numcurvednodes = elementobject.countcurvednodes();
                std::vector<int> nodesincurrentelement(numcurvednodes);
                    
                for (int e = 0; e < elementTags[t].size(); e++)
                {
                    for (int n = 0; n < numcurvednodes; n++)
                    {
                        int currentnode = elemnodeTags[t][e*numcurvednodes+n];
                        if (noderenumbering[currentnode] == -1)
                        {
                            noderenumbering[currentnode] = numberofnodes;
                            numberofnodes++;
                        }
                        nodesincurrentelement[n] = noderenumbering[currentnode];
                    }
                
                    int elementindexincurrenttype = myelements.add(currentelementtype, curvatureorder, nodesincurrentelement);
                    currentphysicalregion->addelement(currentelementtype, elementindexincurrenttype);
                }
            }
        }
    }
    
    
    ///// Create the nodes object:
    
    mynodes.setnumber(numberofnodes);
    std::vector<double>* nodecoordinates = mynodes.getcoordinates();
    
    for (int i = 0; i < noderenumbering.size(); i++)
    {
        int nr = noderenumbering[i];
        int ntr = nodeTagsrenumbering[i];

        if (nr != -1)
        {
            nodecoordinates->at(3*nr+0) = coords[3*ntr+0];
            nodecoordinates->at(3*nr+1) = coords[3*ntr+1];
            nodecoordinates->at(3*nr+2) = coords[3*ntr+2];
        }
    }
}

void gmshinterface::readwithapi(std::string name, nodes& mynodes, elements& myelements, physicalregions& myphysicalregions)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 3);
    gmsh::open(name);
    if (gentools::getfileextension(name) == ".geo")
        gmsh::model::mesh::generate(gmsh::model::getDimension());
    readfromapi(mynodes, myelements, myphysicalregions);
    gmsh::finalize();
}
#endif



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
            gentools::osclean(currentline);
            if (currentline == "$MeshFormat")
            {
                std::getline(meshfile, currentline);
                gentools::osclean(currentline);
                // Get version number:
                mystring stringobject(currentline);
                formatversion = std::stod(stringobject.getstringtonextwhitespace());
                break;
            }
        }
        // Give an error if version is not supported:
        if (formatversion >= 4)
        {
            logs log;
            log.msg() << "Error in 'gmshinterface': GMSH format " << formatversion << " is not supported in the native mesh reader." << std::endl;
            log.msg() << "Use the GMSH API, the petsc mesh reader or export as GMSH 2 format." << std::endl;
            log.error();
        }

        // Move to the node section and read the number of nodes:
        int numberofnodes;
        while (std::getline(meshfile, currentline))
        {
            gentools::osclean(currentline);
            if (currentline == "$Nodes")
            {
                std::getline(meshfile, currentline);
                gentools::osclean(currentline);
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
            gentools::osclean(currentline);
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
            gentools::osclean(currentline);
            if (currentline == "$Elements")
            {
                std::getline(meshfile, currentline);
                gentools::osclean(currentline);
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
            gentools::osclean(currentline);
            mystring stringobject(currentline);
            
            // The first number in the line is an integer (the element number). We skip it:
            stringobject.jumptonextwhitespace();
            // Read now the element type number and define the element object:
            int currentcurvedelementtype = convertgmshelementtypenumber(std::stoi(stringobject.getstringtonextwhitespace()));
            element elementobject(currentcurvedelementtype);
            // Get the uncurved element type number:
            int currentelementtype = elementobject.gettypenumber();
            int curvatureorder = elementobject.getcurvatureorder();
            int elemdim = elementobject.getelementdimension();
            
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
            physicalregion* currentphysicalregion = myphysicalregions.get(universe::physregshift*(elemdim+1) + currentphysicalregionnumber);
            // Read now the node number list in the element. The number of nodes depends on the element:
            std::vector<int> nodesincurrentelement(elementobject.countcurvednodes());
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
        logs log;
        log.msg() << "Unable to open file " << name << " or file not found" << std::endl;
        log.error();
    }
}

void gmshinterface::writetofile(std::string name, nodes& mynodes, elements& myelements, physicalregions& myphysicalregions, disjointregions& mydisjointregions, std::vector<int> physicalregionstowrite)
{    
    // Check which nodes should be written:
    std::vector< std::vector<std::vector<int>>* > elementlists(physicalregionstowrite.size());
    for (int p = 0; p < physicalregionstowrite.size(); p++)
        elementlists[p] = myphysicalregions.get(physicalregionstowrite[p])->getelementlist();
    std::vector<bool> isnodeinphysicalregions;
    int numnodesinphysregs = myelements.istypeinelementlists(0, elementlists, isnodeinphysicalregions, true);
    
    if (numnodesinphysregs == 0)
        return;
        
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
        meshfile << numnodesinphysregs << "\n";
        // Write the node coordinates:        
        std::vector<double>* nodecoordinates = mynodes.getcoordinates();
        
        int index = 0;
        std::vector<int> noderenumbering(mynodes.count(), -1);
        for (int i = 0; i < mynodes.count(); i++)
        {
            if (isnodeinphysicalregions[i])
            {
                meshfile << index+1 << " " << nodecoordinates->at(3*i+0) << " " << nodecoordinates->at(3*i+1) << " " << nodecoordinates->at(3*i+2) << "\n";
                
                noderenumbering[i] = index;
                index++;
            }
        }
        meshfile << "$EndNodes\n";
        
        // Write the element section:
        meshfile << "$Elements\n";
        
        // Write the number of elements.
        // The number of elements is equal to the total number of elements in the physical regions.
        // Elements that are in several physical regions at the same time are counted 
        // multiple times, which is ok in the .msh format.
        int numberofelements = 0;
        for (int p = 0; p < physicalregionstowrite.size(); p++)
            numberofelements += myphysicalregions.get(physicalregionstowrite[p])->countelements();
        meshfile << numberofelements << "\n";
        
        // Write all element lines to the file.
        int elementnumberinfile = 1;
        // Iterate through all physical regions:
        for (int p = 0; p < physicalregionstowrite.size(); p++)
        {
            // Get the physical region object:
            physicalregion* physicalregionobject = myphysicalregions.get(physicalregionstowrite[p]);
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
                    meshfile << elementnumberinfile << " " << converttogmshelementtypenumber(curvedtypenumber) << " 2 " << physicalregionstowrite[p] << " " << physicalregionstowrite[p];
                    // Write all nodes in the element.
                    for (int nodeindex = 0; nodeindex < numberofcurvednodes; nodeindex++)
                        meshfile << " " << noderenumbering[myelements.getsubelement(0, typenumber, i, nodeindex)] + 1; // +1 to start numbering nodes at 1
                    
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
        logs log;
        log.msg() << "Unable to write to file " << name << " or file not found" << std::endl;
        log.error();
    }
}

void gmshinterface::writetofile(std::string name, iodata datatowrite)
{
    // Get the file name without the path and the .pos extension:
    std::string viewname = gentools::getfilename(name);
    
    // Get the list of element types in the view:
    std::vector<int> activeelementtypes = datatowrite.getactiveelementtypes();
    
    // Loop on all active element types:
    for (int i = 0; i < activeelementtypes.size(); i++)
    {
        int elemtypenum = activeelementtypes[i];
        element myelement(elemtypenum);
        
        std::vector<densemat> curcoords = datatowrite.getcoordinates(elemtypenum);
        std::vector<densemat> curdata = datatowrite.getdata(elemtypenum);
        
        // Open the view (overwrite if first time):
        if (activeelementtypes.size() == 1)
            openview(name, viewname, 0, i == 0);
        else
            openview(name, viewname + myelement.gettypename(), 0, i == 0);    
        // Append the data to the view:
        if (datatowrite.isscalar())
            appendtoview(name, elemtypenum, curcoords[0], curcoords[1], curcoords[2], curdata[0]);
        else
            appendtoview(name, elemtypenum, curcoords[0], curcoords[1], curcoords[2], curdata[0], curdata[1], curdata[2]);
        
        // Write the shape function polynomials:
        lagrangeformfunction mylagrange(elemtypenum, datatowrite.getinterpolorder(), {});
        std::vector<polynomial> poly = mylagrange.getformfunctionpolynomials();
        
        std::vector<polynomial> polygeo;
        if (datatowrite.getinterpolorder() == datatowrite.getgeointerpolorder())
            polygeo = poly;
        else
        {
            lagrangeformfunction mylagrangegeo(elemtypenum, datatowrite.getgeointerpolorder(), {});
            polygeo = mylagrangegeo.getformfunctionpolynomials();
        }
        
        writeinterpolationscheme(name, {poly, polygeo});
        // Close the view:
        closeview(name);
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
        logs log;
        log.msg() << "Unable to write to file " << name << " or file not found" << std::endl;
        log.error();
    }
}

void gmshinterface::appendtoview(std::string name, int elementtypenumber, densemat coordx, densemat coordy, densemat coordz, densemat compxinterpolated)
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
        logs log;
        log.msg() << "Unable to write to file " << name << " or file not found" << std::endl;
        log.error();
    }
}

void gmshinterface::appendtoview(std::string name, int elementtypenumber, densemat coordx, densemat coordy, densemat coordz, densemat compxinterpolated, densemat compyinterpolated, densemat compzinterpolated)
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
        logs log;
        log.msg() << "Unable to write to file " << name << " or file not found" << std::endl;
        log.error();
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
        logs log;
        log.msg() << "Unable to write to file " << name << " or file not found" << std::endl;
        log.error();
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
        logs log;
        log.msg() << "Unable to write to file " << name << " or file not found" << std::endl;
        log.error();
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
                logs log;
                log.msg() << "Error in 'gmshinterface' namespace: trying to use a GMSH element (" << gmshtypenumber << ") that is undefined in this code." << std::endl;
                log.error();
        }
    }
    
    throw std::runtime_error(""); // fix return warning
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
                logs log;
                log.msg() << "Error in 'gmshinterface' namespace: trying to use a GMSH element that is undefined in this code." << std::endl;
                log.error();
        }
    }
    
    throw std::runtime_error(""); // fix return warning
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
    
    throw std::runtime_error(""); // fix return warning
}

