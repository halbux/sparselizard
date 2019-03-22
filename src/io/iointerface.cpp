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
    }
    
    std::cout << "Error in 'iointerface' namespace: cannot write to file '" << filename << "' (unknown or missing file extension)" << std::endl;
    abort();
}

bool iointerface::isonlyisoparametric(std::string filename)
{
	// Get the extension:
	std::string fileext = filename.substr(filename.size()-4,4);
	
	if (fileext == ".pos")
		return false;
	if (fileext == ".vtk")
		return true;
		
    std::cout << "Error in 'iointerface' namespace: unknown extension in file '" << filename << std::endl;
    abort();
}
