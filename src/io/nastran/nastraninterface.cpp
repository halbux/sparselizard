#include "nastraninterface.h"
#include "universe.h"


void nastraninterface::readfromfile(std::string name, nodes& mynodes, elements& myelements, physicalregions& myphysicalregions)
{    
    std::string currentline;

    // 'file' cannot take a std::string argument --> name.c_str():
    std::ifstream meshfile (name.c_str());
    if (meshfile.is_open())
    {
        std::vector<std::vector<double>> nodecoords = {};
        
        // Read all GRID and ELEMENT data lines:
        nasdataline mydataline;
        while (std::getline(meshfile, currentline))
        {
            gentools::osclean(currentline);
            bool isready = mydataline.addline(currentline);

            if (isready == true)
            {                
                if (mydataline.isgriddata())
                    nodecoords.push_back(mydataline.getnodecoordinates());
                
                if (mydataline.iselementdata())
                {
                    std::vector<int> nodesincurrentelement = mydataline.getvertices();
                    int curvedelemtypenum = mydataline.getelementtypenumber();
                    int curphysregnum = mydataline.getgroupnumber();
                    
                    element myelem(curvedelemtypenum);
                    int curelemtypenum = myelem.gettypenumber();
                    int elemdim = myelem.getelementdimension();

                    physicalregion* currentphysicalregion = myphysicalregions.get(universe::physregshift*(elemdim+1) + curphysregnum);

                    // Add the element and its physical region:
                    int elementindexincurrenttype = myelements.add(curelemtypenum, myelem.getcurvatureorder(), nodesincurrentelement);
                    currentphysicalregion->addelement(curelemtypenum, elementindexincurrenttype);
                }
            }
        }
        
        
        // Populate the nodes object:
        mynodes.setnumber(nodecoords.size());
        std::vector<double>* nodecoordinates = mynodes.getcoordinates();
        
        for (int i = 0; i < nodecoords.size(); i++)
        {
            nodecoordinates->at(3*i+0) = nodecoords[i][0];
            nodecoordinates->at(3*i+1) = nodecoords[i][1];
            nodecoordinates->at(3*i+2) = nodecoords[i][2];
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

