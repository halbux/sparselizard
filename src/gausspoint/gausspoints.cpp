#include "gausspoints.h"


gausspoints::gausspoints(int elementtypenumber, int integrationorder)
{
    myelementtypenumber = elementtypenumber;
    
    if (integrationorder < 0)
    {
        std::cout << "Error in 'gausspoints' object: cannot get the Gauss points for negative integration order " << integrationorder << std::endl;
        abort();
    }
    
    switch (elementtypenumber)
    {
        // Point element:
        case 0:
            gppoint::set(integrationorder, mycoordinates, myweights);
            break;
        // Line element:
        case 1:
            gpline::set(integrationorder, mycoordinates, myweights);
            break; 
        // Triangle element:
        case 2:
            gptriangle::set(integrationorder, mycoordinates, myweights);
            break; 
        // Quadrangle element:
        case 3:
            gpquadrangle::set(integrationorder, mycoordinates, myweights);
            break; 
        // Tetrahedron element:
        case 4:
            gptetrahedron::set(integrationorder, mycoordinates, myweights);
            break; 
        // Hexahedron element:
        case 5:
            gphexahedron::set(integrationorder, mycoordinates, myweights);
            break; 
        // Prism element:
        case 6:
            gpprism::set(integrationorder, mycoordinates, myweights);
            break; 
        // Pyramid element:
        case 7:
            gppyramid::set(integrationorder, mycoordinates, myweights);
            break; 
        default:
            std::cout << "Error in 'gausspoints' object: unknown element type number " << elementtypenumber << std::endl;
            abort();
    }
}

gausspoints::gausspoints(int elementtypenumber, std::vector<double>& gpcoords)
{
    double roundoffnoise = 1e-12;
    
    myelementtypenumber = elementtypenumber;
    
    if (myelementtypenumber == 0)
    {
        if (gpcoords.size() != 3)
        {
            std::cout << "Error in 'gausspoints' object: could not match gauss coordinates of element type number 0 to the requested coordinates" << std::endl;
            for (int i = 0; i < gpcoords.size(); i++)
                std::cout << gpcoords[i] << " ";
            std::cout << std::endl;
            abort();
        }
    
        gppoint::set(0, mycoordinates, myweights);
        return;
    }
    
    // Find the corresponding integration order or throw an error if none:
    int integrationorder = 0;
    while (true)
    {
        int numgps;
        
        switch (elementtypenumber)
        {
            // Line element:
            case 1:
                numgps = gpline::count(integrationorder);
                if (3*numgps == gpcoords.size())
                    gpline::set(integrationorder, mycoordinates, myweights);
                break; 
            // Triangle element:
            case 2:
                numgps = gptriangle::count(integrationorder);
                if (3*numgps == gpcoords.size())
                    gptriangle::set(integrationorder, mycoordinates, myweights);
                break; 
            // Quadrangle element:
            case 3:
                numgps = gpquadrangle::count(integrationorder);
                if (3*numgps == gpcoords.size())
                    gpquadrangle::set(integrationorder, mycoordinates, myweights);
                break; 
            // Tetrahedron element:
            case 4:
                numgps = gptetrahedron::count(integrationorder);
                if (3*numgps == gpcoords.size())
                    gptetrahedron::set(integrationorder, mycoordinates, myweights);
                break; 
            // Hexahedron element:
            case 5:
                numgps = gphexahedron::count(integrationorder);
                if (3*numgps == gpcoords.size())
                    gphexahedron::set(integrationorder, mycoordinates, myweights);
                break; 
            // Prism element:
            case 6:
                numgps = gpprism::count(integrationorder);
                if (3*numgps == gpcoords.size())
                    gpprism::set(integrationorder, mycoordinates, myweights);
                break; 
            // Pyramid element:
            case 7:
                numgps = gppyramid::count(integrationorder);
                if (3*numgps == gpcoords.size())
                    gppyramid::set(integrationorder, mycoordinates, myweights);
                break; 
            default:
                std::cout << "Error in 'gausspoints' object: unknown element type number " << elementtypenumber << std::endl;
                abort();
        }
        
        if (numgps == -1)
        {
            std::cout << "Error in 'gausspoints' object: could not match gauss coordinates of element type number " << elementtypenumber << " to the requested coordinates" << std::endl;
            for (int i = 0; i < gpcoords.size(); i++)
                std::cout << gpcoords[i] << " ";
            std::cout << std::endl;
            abort();
        }
        
        // Check if close enough to argument coordinates:
        if (3*numgps == gpcoords.size())
        {
            bool iscloseenough = true;
            for (int i = 0; i < gpcoords.size(); i++)
            {
                if (std::abs(gpcoords[i] - mycoordinates[i]) > roundoffnoise)
                {
                    iscloseenough = false;
                    break;
                }
            }
            
            if (iscloseenough)
                break;
        }
        
        integrationorder++;
    }
}

void gausspoints::print(void)
{
    std::cout << std::endl << "ki coordinate | eta coordinate | phi coordinate | weight" << std::endl << std::endl;
    int oldprecision = std::cout.precision();
    std::cout.precision(17);
    for (int i = 0; i < count(); i++)
        std::cout << std::setw(26) << std::left << mycoordinates[3*i+0] << std::setw(26) << std::left << mycoordinates[3*i+1] << std::setw(26) << std::left << mycoordinates[3*i+2] << std::setw(26) << std::left << myweights[i] << std::endl;
    std::cout << std::endl;
    
    std::cout.precision(oldprecision);
}

