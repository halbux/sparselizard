#include "iointerface.h"


void iointerface::writetofile(std::string filename, iodata datatowrite, std::string appendtofilename)
{
    if (filename.size() >= 5)
    {
        // Get the extension:
        std::string fileext = filename.substr(filename.size()-4,4);
        // Get the file name without the extension:
        std::string noext = filename.substr(0, filename.size()-4);
        // Get the requested filename:
        std::string requestedfilename = noext + appendtofilename + fileext;

        if (fileext == ".pos")
        {
            gmshinterface::writetofile(requestedfilename, datatowrite);
            return;
        }
        if (fileext == ".vtk")
        {
            pvinterface::writetovtkfile(requestedfilename, datatowrite);
            return;
        }
        if (fileext == ".vtu")
        {
            pvinterface::writetovtufile(requestedfilename, datatowrite);
            return;
        }
    }
    
    std::cout << "Error in 'iointerface' namespace: cannot write to file '" << filename << "'." << std::endl;
    std::cout << "Supported output formats are .vtk (ParaView), .vtu (ParaView) and .pos (GMSH)." << std::endl;
    abort();
}

bool iointerface::isonlyisoparametric(std::string filename)
{
    if (filename.size() >= 5)
    {
        // Get the extension:
        std::string fileext = filename.substr(filename.size()-4,4);
        
        if (fileext == ".pos")
            return false;
        if (fileext == ".vtk")
            return true;
        if (fileext == ".vtu")
            return true;
    }
		
    std::cout << "Error in 'iointerface' namespace: cannot handle file '" << filename << "'." << std::endl;
    std::cout << "Supported output formats are .vtk (ParaView), .vtu (ParaView) and .pos (GMSH)." << std::endl;
    abort();
}

void iointerface::write(std::string filename, std::vector<int>& intdata, std::vector<double>& doubledata, bool isbinary)
{
    if (isbinary == false)
    {
        // 'file' cannot take a std::string argument --> filename.c_str():
        std::ofstream name (filename.c_str());
        if (name.is_open())
        {
            // Write the int and double data length:
            name << intdata.size() << std::endl << doubledata.size() << std::endl;

            // Write the int data to the file:
            for (int i = 0; i < intdata.size(); i++)
                name << intdata[i] << std::endl;

            // To write all doubles with enough digits to the file:
            name << std::setprecision(17);

            for (int i = 0; i < doubledata.size()-1; i++)
                name << doubledata[i] << std::endl;
            name << doubledata[doubledata.size()-1];

            name.close();
        }
        else 
        {
            std::cout << "Unable to write data to file " << filename << " or file not found" << std::endl;
            abort();
        } 
    }
    else
    {
        // All ints will be converted to doubles (doubles can exactly represent ints up to at least 2^52):
        int totalsize = intdata.size() + doubledata.size() + 2;
        intdensematrix alladdresses(1,totalsize, 0,1);
        densematrix alldata(1,totalsize);

        int* addsvals = alladdresses.getvalues();
        double* datavals = alldata.getvalues();

        datavals[0] = intdata.size();
        datavals[1] = doubledata.size();
        for (int i = 0; i < intdata.size(); i++)
            datavals[2 + i] = intdata[i];
        for (int i = 0; i < doubledata.size(); i++)
            datavals[2+intdata.size() + i] = doubledata[i];
            
        // Let petsc write the binary file for us:
        Vec datvec;
        VecCreate(PETSC_COMM_SELF, &datvec);
        VecSetSizes(datvec, PETSC_DECIDE, totalsize);
        VecSetFromOptions(datvec);  

        VecSetValues(datvec, totalsize, addsvals, datavals, INSERT_VALUES);

        PetscViewer v;
        PetscViewerBinaryOpen(PETSC_COMM_SELF, filename.c_str(), FILE_MODE_WRITE, &v);
        VecView(datvec, v);
        PetscViewerDestroy(&v);

        VecDestroy(&datvec);
    }
}

void iointerface::load(std::string filename, std::vector<int>& intdata, std::vector<double>& doubledata, bool isbinary)
{
    if (isbinary == false)
    {
        std::string currentline;

        // 'file' cannot take a std::string argument --> filename.c_str():
        std::ifstream name (filename.c_str());
        if (name.is_open())
        {
            // First integer is the intdata size, second is the doubledata size:
            std::getline(name, currentline);
            intdata.resize(std::stoi(currentline));
            
            std::getline(name, currentline);
            doubledata.resize(std::stoi(currentline));

            // Load the intdata:
            for (int i = 0; i < intdata.size(); i++)
            {
                std::getline(name, currentline);
                intdata[i] = std::stoi(currentline);
            }
            
            // Load the doubledata:
            for (int i = 0; i < doubledata.size(); i++)
            {
                std::getline(name, currentline);
                doubledata[i] = std::stod(currentline);
            }

            name.close();
        }
        else 
        {
            std::cout << "Unable to load data from file " << filename << " or file not found" << std::endl;
            abort();
        }
    }
    else
    {
        // Let petsc read the binary file for us:
        Vec datvec;
        VecCreate(PETSC_COMM_SELF, &datvec);
        
        PetscViewer v;
        PetscViewerBinaryOpen(PETSC_COMM_SELF, filename.c_str(), FILE_MODE_READ, &v);
        VecLoad(datvec, v);
        PetscViewerDestroy(&v);
        
        int veclen;
        VecGetSize(datvec, &veclen);
        
        densematrix doublestoget(1, veclen);
        intdensematrix addressestoget(1, veclen, 0, 1);
        
        double* vals = doublestoget.getvalues();
        int* ads = addressestoget.getvalues();
        
        VecGetValues(datvec, veclen, ads, vals);
        
        int numints = vals[0];
        int numdoubles = vals[1];
        intdata.resize(numints);
        doubledata.resize(numdoubles);
        
        for (int i = 0; i < numints; i++)
            intdata[i] = vals[2 + i];
        for (int i = 0; i < numdoubles; i++)
            doubledata[i] = vals[2+numints + i];
        
        VecDestroy(&datvec);
    }
}

