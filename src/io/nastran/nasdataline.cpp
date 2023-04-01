#include "nasdataline.h"


std::vector<int> nasdataline::translateelementname(std::string elemname)
{
    if (elemname == "CBAR" || elemname == "CROD" || elemname == "CBEAM")
        return {1,2};
    if (elemname == "CTRIA3")
        return {2,3};
    if (elemname == "CQUAD4")
        return {3,4};
    if (elemname == "CTETRA")
        return {4,4};
    if (elemname == "CHEXA")
        return {5,8};
    if (elemname == "CPENTA")
        return {6,6};
        
    logs log;
    log.msg() << "Error in 'nasdataline': unknown or unsupported Nastran element type '" << elemname << "'." << std::endl;
    log.msg() << "Curved elements not supported in .nas reader (save as linear elements in .nas or use GMSH 2 ASCII .msh format for curved elements)." << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

bool nasdataline::addline(std::string linetoadd)
{
    // In case we are at a new 'GRID' line:
    if (linetoadd.length() >= 4 && linetoadd.compare(0,4,"GRID") == 0)
    {
        isitgriddata = true; isitelementdata = false;

        // Check for format:
        if (linetoadd[4] == ' ')
            lineformat = 1;
        if (linetoadd[4] == '*')
            lineformat = 2;
        if (linetoadd[4] == ',')
            lineformat = 3;

        if (lineformat == 1)
        {
            nodecoordinate[0] = std::stod(linetoadd.substr(24,8));
            nodecoordinate[1] = std::stod(linetoadd.substr(32,8));
            nodecoordinate[2] = std::stod(linetoadd.substr(40,8));
            return true; 
        }
        if (lineformat == 2)
        {
            nodecoordinate[0] = std::stod(linetoadd.substr(40,16));
            nodecoordinate[1] = std::stod(linetoadd.substr(56,16));
            // Z coordinate does not fit in this line, will be on next one.
        }
        if (lineformat == 3)
        {
            mystring curline(linetoadd);
            curline.getstringtonextcomma();
            curline.getstringtonextcomma();
            curline.getstringtonextcomma();
            nodecoordinate[0] = std::stod(curline.getstringtonextcomma());
            nodecoordinate[1] = std::stod(curline.getstringtonextcomma());
            nodecoordinate[2] = std::stod(curline.getstringtonextcomma());
            return true;
        }
    }

    // In case we are continuing the previous 'GRID' line:
    if (isitgriddata && linetoadd.length() > 1 && linetoadd[0] == '*')
    {
        nodecoordinate[2] = std::stod(linetoadd.substr(8,16));
        return true;
    }
    
    
    // In case we are at a new element data line:
    if (linetoadd.length() >= 2 && linetoadd.compare(0,1,"C") == 0)
    {
        isitgriddata = false; isitelementdata = true;

        mystring curline(linetoadd);
            
        if (lineformat == 1 || lineformat == 2)
        {
            std::vector<int> eleminfo = translateelementname(curline.getstringtonextwhitespace());
            elementtypenumber = eleminfo[0];
            groupnumber = std::stoi(linetoadd.substr(16,8));
            
            vertices.resize(eleminfo[1]);
            
            for (int i = 0; i < vertices.size(); i++)
            {
                std::string curstrng = linetoadd.substr(24+8*i,8);
                if (curstrng[0] == '+')
                    break;
                vertices[i] = std::stoi(curstrng) - 1;
                currentvertexindex = i+1;
            }
            if (currentvertexindex == vertices.size())
            {
                currentvertexindex = 0;
                return true;
            }
            else
                return false;
        }
        if (lineformat == 3)
        {
            std::vector<int> eleminfo = translateelementname(curline.getstringtonextcomma());
            elementtypenumber = eleminfo[0];
            curline.getstringtonextcomma();
            groupnumber = std::stoi(curline.getstringtonextcomma());
            
            vertices.resize(eleminfo[1]);

            for (int i = 0; i < vertices.size(); i++)
            {
                std::string curstrng = curline.getstringtonextcomma();
                if (curstrng[0] == '+')
                    break;
                vertices[i] = std::stoi(curstrng) - 1;
                currentvertexindex = i+1;
            }
            if (currentvertexindex == vertices.size())
            {
                currentvertexindex = 0;
                return true;
            }
            else
                return false;
        }
    }

    // In case we are continuing the previous element data line:
    if (isitelementdata && linetoadd.length() > 1 && linetoadd[0] == '+')
    {
        mystring curline(linetoadd);
    
        if (lineformat == 1 || lineformat == 2)
        {
            curline.getstringtonextwhitespace();

            int i;
            for (i = 0; i < vertices.size()-currentvertexindex; i++)
            {
                std::string curstrng = linetoadd.substr(8+8*i,8);
                if (curstrng[0] == '+')
                    break;
                vertices[currentvertexindex+i] = std::stoi(curstrng) - 1;
            }
            currentvertexindex += i;
            if (currentvertexindex == vertices.size())
            {
                currentvertexindex = 0;
                return true;
            }
            else
                return false;
        }
        if (lineformat == 3)
        {
            curline.getstringtonextcomma();

            int i;
            for (i = 0; i < vertices.size()-currentvertexindex; i++)
            {
                std::string curstrng = curline.getstringtonextcomma();
                if (curstrng[0] == '+')
                    break;
                vertices[currentvertexindex+i] = std::stoi(curstrng) - 1;
            }
            currentvertexindex += i;
            if (currentvertexindex == vertices.size())
            {
                currentvertexindex = 0;
                return true;
            }
            else
                return false;
        }
    }


    return false;
}

