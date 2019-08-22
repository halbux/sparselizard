#include "gausspoints.h"


gausspoints::gausspoints(int elementtypenumber, int integrationorder)
{
    myelementtypenumber = elementtypenumber;
    myintegrationorder = integrationorder;
    
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

