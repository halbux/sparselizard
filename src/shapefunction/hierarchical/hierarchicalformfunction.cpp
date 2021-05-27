#include "hierarchicalformfunction.h"


int hierarchicalformfunction::gettypenumber(std::string fftypename)
{
    int outtypenum = -1;
    for (int i = 0; i < mytypenames.size(); i++)
    {
        if (mytypenames[i] == fftypename)
        {
            outtypenum = i;
            break;
        }
    }
    
    if (outtypenum == -1)
    {
        std::cout << "Error in 'hierarchicalformfunction' object: unknown type name '" << fftypename << "'" << std::endl;
        abort();
    }
    
    return outtypenum;
}

std::string hierarchicalformfunction::gettypename(int fftypenumber)
{
    if (fftypenumber >= mytypenames.size())
    {
        std::cout << "Error in 'hierarchicalformfunction' object: unknown type number " << fftypenumber << std::endl;
        abort();
    }

    return mytypenames[fftypenumber];
}

int hierarchicalformfunction::getminorder(std::string fftypename)
{
    if (fftypename == "hcurl")
        return 0;
    if (fftypename == "h1" || fftypename == "one0" || fftypename == "one1" || fftypename == "one2" || fftypename == "one3" || fftypename == "h1d0" || fftypename == "h1d1" || fftypename == "h1d2" || fftypename == "h1d3" || fftypename == "x" || fftypename == "y" || fftypename == "z")
        return 1;
        
    std::cout << "Error in 'hierarchicalformfunction' object: unknown type name '" << fftypename << "'" << std::endl;
    abort();
}

int hierarchicalformfunction::getminorderdim(std::string fftypename)
{
    if (fftypename == "h1" || fftypename == "one0" || fftypename == "h1d0" || fftypename == "x" || fftypename == "y" || fftypename == "z")
        return 0;
    if (fftypename == "hcurl" || fftypename == "one1" || fftypename == "h1d1")
        return 1;
    if (fftypename == "one2" || fftypename == "h1d2")
        return 2;
    if (fftypename == "one3" || fftypename == "h1d3")
        return 3;
        
    std::cout << "Error in 'hierarchicalformfunction' object: unknown type name '" << fftypename << "'" << std::endl;
    abort();
}

