#include "iointerface.h"


void iointerface::writetofile(std::string filename, iodata datatowrite)
{
    if (filename.size() >= 5)
    {
        // Get the extension:
        std::string fileext = filename.substr(filename.size()-4,4);
        
        if (fileext == ".pos")
        {
            gmshinterface::writetofile(filename, datatowrite);
            return;
        }
    }
    
    std::cout << "Error in 'iointerface' namespace: cannot write to file '" << filename << "' (unknown or missing file extension)" << std::endl;
    abort();
}
