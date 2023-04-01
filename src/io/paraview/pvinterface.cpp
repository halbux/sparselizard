#include "pvinterface.h"


void pvinterface::writetovtkfile(std::string name, iodata datatowrite)
{
    // Get the file name without the .vtk extension:
    std::string namenoext = name.substr(0, name.size()-4);

    // Get all timesteps in 'datatowrite':
    std::vector<double> timetags = datatowrite.gettimetags();

    if (timetags.size() == 0)
        writetovtkfile(name, datatowrite, -1);

    for (int i = 0; i < timetags.size(); i++)
    {
        std::string curname = namenoext + "_" + std::to_string(i) + ".vtk";
        writetovtkfile(curname, datatowrite, i);
    }
}

void pvinterface::writetovtufile(std::string name, iodata datatowrite)
{
    // Get the file name without the .vtu extension:
    std::string namenoext = name.substr(0, name.size()-4);

    // Get all timesteps in 'datatowrite':
    std::vector<double> timetags = datatowrite.gettimetags();

    if (timetags.size() == 0)
        writetovtufile(name, datatowrite, -1);

    for (int i = 0; i < timetags.size(); i++)
    {
        std::string curname = namenoext + "_" + std::to_string(i) + ".vtu";
        writetovtufile(curname, datatowrite, i);
    }
}

void pvinterface::writetovtkfile(std::string name, iodata datatowrite, int timestepindex)
{
    // Get the file name without the path and the .vtk extension:
    std::string viewname = gentools::getfilename(name);

    mystring myname(viewname);
    viewname = myname.getstringwhileletter();

    // 'file' cannot take a std::string argument --> name.c_str():
    std::ofstream outfile (name.c_str());
    if (outfile.is_open())
    {
        // To write all doubles with enough digits to the file:
        outfile << std::setprecision(17);

        // Write the header:
        outfile << "# vtk DataFile Version 4.2\n";
        outfile << viewname+"\n";
        outfile << "ASCII\n";
        outfile << "DATASET UNSTRUCTURED_GRID\n\n";

        // Write the points section.
        int numnodes = datatowrite.countcoordnodes();
        outfile << "POINTS " << numnodes << " double\n";
        for (int tn = 0; tn < 8; tn++)
        {
            if (datatowrite.ispopulated(tn) == false)
                continue;

            std::vector<densemat> curcoords = datatowrite.getcoordinates(tn,timestepindex);
            double* xvals = curcoords[0].getvalues();
            double* yvals = curcoords[1].getvalues();
            double* zvals = curcoords[2].getvalues();

            int index = 0;
            for (int elem = 0; elem < curcoords[0].countrows(); elem++)
            {
                for (int node = 0; node < curcoords[0].countcolumns(); node++)
                {
                    if (universe::isaxisymmetric)
                        outfile << xvals[index] << " " << -zvals[index] << " " << yvals[index] << " ";
                    else
                        outfile << xvals[index] << " " << yvals[index] << " " << zvals[index] << " ";
                
                    index++;
                }
                outfile << "\n";
            }
        }
        outfile << "\n";

        // Write the cells section:
        int numelems = datatowrite.countelements();
        outfile << "CELLS " << numelems << " " << numnodes + numelems << "\n";

        int nodenum = 0;
        for (int tn = 0; tn < 8; tn++)
        {
            if (datatowrite.ispopulated(tn) == false)
                continue;

            // Move from our node ordering to the one of ParaView:
            element myelem(tn, datatowrite.getinterpolorder());
            std::vector<int> reordering = getnodereordering(myelem.getcurvedtypenumber());

            std::vector<densemat> curcoords = datatowrite.getcoordinates(tn,timestepindex);

            for (int elem = 0; elem < curcoords[0].countrows(); elem++)
            {
                outfile << curcoords[0].countcolumns() << " ";
                for (int node = 0; node < curcoords[0].countcolumns(); node++)
                    outfile << nodenum + reordering[node] << " ";
                nodenum += curcoords[0].countcolumns();
                outfile << "\n";
            }
        }
        outfile << "\n";

        // Write the cell types section:
        outfile << "CELL_TYPES " << numelems << "\n";
        for (int tn = 0; tn < 8; tn++)
        {
            if (datatowrite.ispopulated(tn) == false)
                continue;

            element myelem(tn, datatowrite.getinterpolorder());

            std::vector<densemat> curcoords = datatowrite.getcoordinates(tn,timestepindex);

            for (int elem = 0; elem < curcoords[0].countrows(); elem++)
                outfile << converttoparaviewelementtypenumber(myelem.getcurvedtypenumber()) << "\n";
        }
        outfile << "\n";

        // Write the data section:
        outfile << "POINT_DATA " << numnodes << "\n";
        outfile << "\n";

        // Write the scalar data section (if any):
        if (datatowrite.isscalar() == true)
        {
            outfile << "SCALARS " << viewname << " double\n";
            outfile << "LOOKUP_TABLE default" << "\n";
            for (int tn = 0; tn < 8; tn++)
            {
                if (datatowrite.ispopulated(tn) == false)
                    continue;

                densemat scaldat = datatowrite.getdata(tn,timestepindex)[0];
                double* scalvals = scaldat.getvalues();

                for (int i = 0; i < scaldat.count(); i++)
                    outfile << scalvals[i] << "\n";
            }
            outfile << "\n";
        }

        // Write the vector data section (if any):
        if (datatowrite.isscalar() == false)
        {
            outfile << "VECTORS " << viewname << " double\n";
            for (int tn = 0; tn < 8; tn++)
            {
                if (datatowrite.ispopulated(tn) == false)
                    continue;

                std::vector<densemat> vecdat = datatowrite.getdata(tn,timestepindex);
                double* compxvals = vecdat[0].getvalues();
                double* compyvals = vecdat[1].getvalues();
                double* compzvals = vecdat[2].getvalues();

                for (int i = 0; i < vecdat[0].count(); i++)
                {
                    if (universe::isaxisymmetric)
                        outfile << compxvals[i] << " " << -compzvals[i]  << " " << compyvals[i] << "\n";
                    else
                        outfile << compxvals[i] << " " << compyvals[i] << " " << compzvals[i] << "\n";
                }
            }
            outfile << "\n";
        }

        outfile.close();
    }
    else 
    {
        logs log;
        log.msg() << "Unable to write to file " << name << " or file not found" << std::endl;
        log.error();
    }
}

void pvinterface::writetovtufile(std::string name, iodata datatowrite, int timestepindex)
{
    // Get the file name without the path and the .vtu extension:
    std::string viewname = gentools::getfilename(name);

    mystring myname(viewname);
    viewname = myname.getstringwhileletter();

    // 'file' cannot take a std::string argument --> name.c_str():
    std::ofstream outfile (name.c_str());
    if (outfile.is_open())
    {
        // To write all doubles with enough digits to the file:
        outfile << std::setprecision(17);
        
        // Write the header:
        outfile << "<?xml version=\"1.0\"?>\n";
        outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        outfile << "<UnstructuredGrid>\n";

        // Get the number of points and cells
        int numnodes = datatowrite.countcoordnodes();
        int numelems = datatowrite.countelements();
        outfile << "<Piece NumberOfPoints=\"" << numnodes << "\" NumberOfCells=\"" << numelems << "\">\n";
            
        // Write the points section.
        outfile << "<Points>\n";
        outfile << "<DataArray type=\"Float64\" Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        
        for (int tn = 0; tn < 8; tn++)
        {
            if (datatowrite.ispopulated(tn) == false)
                continue;

            std::vector<densemat> curcoords = datatowrite.getcoordinates(tn,timestepindex);
            double* xvals = curcoords[0].getvalues();
            double* yvals = curcoords[1].getvalues();
            double* zvals = curcoords[2].getvalues();

            int index = 0;
            for (int elem = 0; elem < curcoords[0].countrows(); elem++)
            {
                for (int node = 0; node < curcoords[0].countcolumns(); node++)
                {
                    if (universe::isaxisymmetric)
                        outfile << xvals[index] << " " << -zvals[index] << " " << yvals[index] << "\n";
                    else
                        outfile << xvals[index] << " " << yvals[index] << " " << zvals[index] << "\n";
                        
                    index++;
                }
            }
        }
        outfile << "</DataArray>\n";
        outfile << "</Points>\n";

        // Write the cells section:
        outfile << "<Cells>\n";
        outfile << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";

        int nodenum = 0;
        for (int tn = 0; tn < 8; tn++)
        {
            if (datatowrite.ispopulated(tn) == false)
                continue;

            // Move from our node ordering to the one of ParaView:
            element myelem(tn, datatowrite.getinterpolorder());
            std::vector<int> reordering = getnodereordering(myelem.getcurvedtypenumber());

            std::vector<densemat> curcoords = datatowrite.getcoordinates(tn,timestepindex);

            for (int elem = 0; elem < curcoords[0].countrows(); elem++)
            {
                for (int node = 0; node < curcoords[0].countcolumns(); node++)
                    outfile << nodenum + reordering[node] << " ";
                nodenum += curcoords[0].countcolumns();
                outfile << "\n";
            }
        }
        outfile << "</DataArray>\n";
        
        // Write the offset section:
        outfile << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        
        int nb = 0;
        for (int tn = 0; tn < 8; tn++)
        {
            if (datatowrite.ispopulated(tn) == false)
                continue;

            std::vector<densemat> curcoords = datatowrite.getcoordinates(tn,timestepindex);

            for (int elem = 0; elem < curcoords[0].countrows(); elem++)
            {
                nb += curcoords[0].countcolumns();
                outfile << nb << "\n";
            }
        }
        outfile << "</DataArray>\n";

        // Write the cell types section:
        outfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        
        for (int tn = 0; tn < 8; tn++)
        {
            if (datatowrite.ispopulated(tn) == false)
                continue;

            element myelem(tn, datatowrite.getinterpolorder());

            std::vector<densemat> curcoords = datatowrite.getcoordinates(tn,timestepindex);

            for (int elem = 0; elem < curcoords[0].countrows(); elem++)
                outfile << converttoparaviewelementtypenumber(myelem.getcurvedtypenumber()) << "\n";
        }
        outfile << "</DataArray>\n";
        outfile << "</Cells>\n";

        // Write the data section:

        // Write the scalar data section (if any):
        if (datatowrite.isscalar() == true)
        {
            outfile << "<PointData Scalars=\"" << viewname << "\">\n";
            outfile << "<DataArray type=\"Float64\" Name=\"" << viewname << "\" format=\"ascii\">\n";
            for (int tn = 0; tn < 8; tn++)
            {
                if (datatowrite.ispopulated(tn) == false)
                    continue;

                densemat scaldat = datatowrite.getdata(tn,timestepindex)[0];
                double* scalvals = scaldat.getvalues();

                for (int i = 0; i < scaldat.count(); i++)
                    outfile << scalvals[i] << "\n";
            }
            outfile << "</DataArray>\n";
            outfile << "</PointData>\n";
        }

        // Write the vector data section (if any):
        if (datatowrite.isscalar() == false)
        {
            outfile << "<PointData Vectors=\"" << viewname << "\">\n";
            outfile << "<DataArray type=\"Float64\" Name=\"" << viewname << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (int tn = 0; tn < 8; tn++)
            {
                if (datatowrite.ispopulated(tn) == false)
                    continue;

                std::vector<densemat> vecdat = datatowrite.getdata(tn,timestepindex);
                double* compxvals = vecdat[0].getvalues();
                double* compyvals = vecdat[1].getvalues();
                double* compzvals = vecdat[2].getvalues();

                for (int i = 0; i < vecdat[0].count(); i++)
                {
                    if (universe::isaxisymmetric)
                        outfile << compxvals[i] << " " << -compzvals[i]  << " " << compyvals[i] << "\n";
                    else
                        outfile << compxvals[i] << " " << compyvals[i] << " " << compzvals[i] << "\n";
                }
            }
            outfile << "</DataArray>\n";
            outfile << "</PointData>\n";
        }

        outfile << "</Piece>\n";
        outfile << "</UnstructuredGrid>\n";
        outfile << "</VTKFile>\n";
            
        outfile.close();
    }
    else 
    {
        logs log;
        log.msg() << "Unable to write to file " << name << " or file not found" << std::endl;
        log.error();
    }
}

void pvinterface::grouptopvdfile(std::string filename, std::vector<std::string> filestogroup, std::vector<double> timevals)
{
    int numsteps = timevals.size();

    // 'file' cannot take a std::string argument --> filename.c_str():
    std::ofstream outfile (filename.c_str());
    if (outfile.is_open())
    {
        // To write all doubles with enough digits to the file:
        outfile << std::setprecision(17);
        
        // Write the header:
        outfile << "<?xml version=\"1.0\"?>\n";
        outfile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        outfile << "<Collection>\n";
    
        for (int i = 0; i < numsteps; i++)
        {
            std::ifstream curfile (filestogroup[i].c_str());
            if (curfile.is_open())
            {
                curfile.close();
                outfile << "<DataSet timestep=\"" << timevals[i] << "\" group=\"\" part=\"0\" file=\"" << filestogroup[i] << "\"/>\n";
            }
            else
            {
                logs log;
                log.msg() << "Error in 'pvinterface': could not find file '" << filestogroup[i] << "' during grouping" << std::endl;
                log.error();
            }
        }
    
        outfile << "</Collection>\n";
        outfile << "</VTKFile>\n";
    
        outfile.close();
    }
    else 
    {
        logs log;
        log.msg() << "Unable to write to file " << filename << " or file not found" << std::endl;
        log.error();
    }
}

int pvinterface::converttoparaviewelementtypenumber(int ourtypenumber)
{
    // Point:
    if (ourtypenumber == 0)
        return 1;

    // This is the general Lagrange elements starting from ParaView version 5.5:
    element myelement(ourtypenumber);
    return (67 + myelement.gettypenumber());
}

std::vector<int> pvinterface::getnodereordering(int ourtypenumber)
{
    element myelement(ourtypenumber);
    
    int order = myelement.getcurvatureorder();
    int elemtypenumber = myelement.gettypenumber();
    
    // Point:
    if (elemtypenumber == 0)
        return {0};

    // Line:
    if (elemtypenumber == 1)
    {
        switch (order)
        {
            case 1:
                return {0,1};
            case 2:
                return {0,1,2};
            case 3:
                return {0,1,2,3};
            case 4:
                return {0,1,2,3,4};
            case 5:
                return {0,1,2,3,4,5};
        }
    }

    // Triangle:
    if (elemtypenumber == 2)
    {
        switch (order)
        {
            case 1:
                return {0,1,2};
            case 2:
                return {0,1,2,3,4,5};
            case 3:
                return {0,1,2,3,4,5,6,7,8,9};
            case 4:
                return {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
            case 5:
                return {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
        }
    }

    // Quadrangle:
    if (elemtypenumber == 3)
    {
        switch (order)
        {
            case 1:
                return {0,1,2,3};
            case 2:
                return {0,1,2,3,4,5,6,7,8};
            case 3:
                return {0,1,2,3,4,5,6,7,9,8,11,10,12,13,15,14};
            case 4:
                return {0,1,2,3,4,5,6,7,8,9,12,11,10,15,14,13,16,20,17,23,24,21,19,22,18};
            case 5:
                return {0,1,2,3,4,5,6,7,8,9,10,11,15,14,13,12,19,18,17,16,20,24,25,21,31,32,33,26,30,35,34,27,23,29,28,22};
        }
    }

    // Tetrahedron:
    if (elemtypenumber == 4)
    {
        switch (order)
        {
            case 1:
                return {0,1,2,3};
            case 2:
                return {0,1,2,3,4,5,6,7,9,8};
            case 3:
                return {0,1,2,3,4,5,6,7,8,9,11,10,15,14,13,12,17,19,18,16};
            case 4:
                return {0,1,2,3,4,5,6,7,8,9,10,11,12,15,14,13,21,20,19,18,17,16,25,26,27,33,31,32,28,29,30,22,23,24,34};
            case 5:
                return {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,19,18,17,16,27,26,25,24,23,22,21,20,34,35,36,37,38,39,48,46,47,51,49,50,40,41,42,43,44,45,28,29,30,31,32,33,52,53,54,55};
        }
    }

    // Hexahedron:
    if (elemtypenumber == 5)
    {
        switch (order)
        {
            case 1:
                return {0,1,2,3,4,5,6,7};
            case 2:
                return {0,1,2,3,4,5,6,7,8,11,13,9,16,18,19,17,10,12,15,14,22,23,21,24,20,25,26};
            case 3:
                return {0,1,2,3,4,5,6,7,8,9,14,15,19,18,10,11,24,25,28,29,31,30,26,27,12,13,16,17,22,23,20,21,40,43,41,42,44,45,47,46,36,37,39,38,49,48,50,51,32,35,33,34,52,53,55,54,56,57,59,58,60,61,63,62};
            case 4:
                return {0,1,2,3,4,5,6,7,8,9,10,17,18,19,25,24,23,11,12,13,32,33,34,38,39,40,43,42,41,35,36,37,14,15,16,20,21,22,29,30,31,26,27,28,62,69,65,66,70,68,63,67,64,71,75,72,78,79,76,74,77,73,53,57,54,60,61,58,56,59,55,81,84,80,85,88,87,82,86,83,44,51,47,48,52,50,45,49,46,89,93,90,96,97,94,92,95,91,98,106,99,107,118,109,101,111,100,108,119,110,120,124,121,113,122,112,102,114,103,115,123,116,105,117,104};
            case 5:
                return {0,1,2,3,4,5,6,7,8,9,10,11,20,21,22,23,31,30,29,28,12,13,14,15,40,41,42,43,48,49,50,51,55,54,53,52,44,45,46,47,16,17,18,19,24,25,26,27,36,37,38,39,32,33,34,35,88,99,98,91,92,100,103,97,93,101,102,96,89,94,95,90,104,108,109,105,115,116,117,110,114,119,118,111,107,113,112,106,72,76,77,73,83,84,85,78,82,87,86,79,75,81,80,74,121,125,124,120,126,133,132,131,127,134,135,130,122,128,129,123,56,67,66,59,60,68,71,65,61,69,70,64,57,62,63,58,136,140,141,137,147,148,149,142,146,151,150,143,139,145,144,138,152,160,161,153,162,184,187,166,163,185,186,167,155,171,170,154,164,188,189,168,192,208,209,196,195,211,210,197,174,201,200,172,165,191,190,169,193,212,213,199,194,215,214,198,175,202,203,173,156,176,177,157,178,204,205,180,179,207,206,181,159,183,182,158};
        }
    }

    // Prism:
    if (elemtypenumber == 6)
    {
        switch (order)
        {
            case 1:
                return {0,1,2,3,4,5};
            case 2:
                return {0,1,2,3,4,5,6,9,7,12,14,13,8,10,11,15,17,16};
            case 3:
                return {0,1,2,3,4,5,6,7,12,13,9,8,18,19,22,23,21,20,10,11,14,15,16,17,24,25,26,27,29,28,34,35,37,36,33,30,32,31,38,39};
            case 4:
                return {0,1,2,3,4,5,6,7,8,15,16,17,11,10,9,24,25,26,30,31,32,29,28,27,12,13,14,18,19,20,21,22,23,33,35,34,36,37,38,39,43,40,46,47,44,42,45,41,57,61,58,64,65,62,60,63,59,51,55,48,54,56,52,50,53,49,66,69,72,68,71,74,67,70,73};
            case 5:
                return {0,1,2,3,4,5,6,7,8,9,18,19,20,21,13,12,11,10,30,31,32,33,38,39,40,41,37,36,35,34,14,15,16,17,22,23,24,25,26,27,28,29,42,47,44,45,46,43,48,51,49,53,52,50,54,58,59,55,65,66,67,60,64,69,68,61,57,63,62,56,86,90,91,87,97,98,99,92,96,101,100,93,89,95,94,88,73,80,81,70,79,85,82,74,78,84,83,75,72,77,76,71,102,114,106,122,118,110,104,116,108,124,120,112,105,117,109,125,121,113,103,115,107,123,119,111};
        }
    }
    
    logs log;
    log.msg() << "Error in 'pvinterface' namespace: trying to use a ParaView element that is undefined in this code." << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}


