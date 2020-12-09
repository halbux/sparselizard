#include "selector.h"


std::shared_ptr<hierarchicalformfunction> selector::select(int elementtypenumber, std::string formfunctiontypename)
{
    if (formfunctiontypename.compare("h1") == 0)
    {
        switch (elementtypenumber)
        {
            case 0:
                return std::shared_ptr<hierarchicalformfunction>(new h1point(-1));
            case 1:
                return std::shared_ptr<hierarchicalformfunction>(new h1line(-1));
            case 2:
                return std::shared_ptr<hierarchicalformfunction>(new h1triangle(-1));
            case 3:
                return std::shared_ptr<hierarchicalformfunction>(new h1quadrangle(-1));
            case 4:
                return std::shared_ptr<hierarchicalformfunction>(new h1tetrahedron(-1));
            case 5:
                return std::shared_ptr<hierarchicalformfunction>(new h1hexahedron(-1));
            case 6:
                return std::shared_ptr<hierarchicalformfunction>(new h1prism(-1));
            case 7:
                return std::shared_ptr<hierarchicalformfunction>(new h1pyramid(-1));
        }
    }
    
    if (formfunctiontypename.compare("hcurl") == 0)
    {
        switch (elementtypenumber)
        {
            case 0:
                return std::shared_ptr<hierarchicalformfunction>(new hcurlpoint);
            case 1:
                return std::shared_ptr<hierarchicalformfunction>(new hcurlline);
            case 2:
                return std::shared_ptr<hierarchicalformfunction>(new hcurltriangle);
            case 3:
                return std::shared_ptr<hierarchicalformfunction>(new hcurlquadrangle);
            case 4:
                return std::shared_ptr<hierarchicalformfunction>(new hcurltetrahedron);
            case 5:
                return std::shared_ptr<hierarchicalformfunction>(new hcurlhexahedron);
            case 6:
                return std::shared_ptr<hierarchicalformfunction>(new hcurlprism);
            case 7:
                return std::shared_ptr<hierarchicalformfunction>(new hcurlpyramid);
        }
    }
    
    if (formfunctiontypename.compare("one0") == 0 || formfunctiontypename.compare("one1") == 0 || formfunctiontypename.compare("one2") == 0 || formfunctiontypename.compare("one3") == 0)
    {
        int onedim = std::stoi(formfunctiontypename.substr(3,1));
        
        return std::shared_ptr<hierarchicalformfunction>(new oneconstant(onedim, elementtypenumber));
    }
    
    if (formfunctiontypename.compare("h1d0") == 0 || formfunctiontypename.compare("h1d1") == 0 || formfunctiontypename.compare("h1d2") == 0 || formfunctiontypename.compare("h1d3") == 0)
    {
        int h1ddim = std::stoi(formfunctiontypename.substr(3,1));
        
        switch (elementtypenumber)
        {
            case 0:
                return std::shared_ptr<hierarchicalformfunction>(new h1point(h1ddim));
            case 1:
                return std::shared_ptr<hierarchicalformfunction>(new h1line(h1ddim));
            case 2:
                return std::shared_ptr<hierarchicalformfunction>(new h1triangle(h1ddim));
            case 3:
                return std::shared_ptr<hierarchicalformfunction>(new h1quadrangle(h1ddim));
            case 4:
                return std::shared_ptr<hierarchicalformfunction>(new h1tetrahedron(h1ddim));
            case 5:
                return std::shared_ptr<hierarchicalformfunction>(new h1hexahedron(h1ddim));
            case 6:
                return std::shared_ptr<hierarchicalformfunction>(new h1prism(h1ddim));
            case 7:
                return std::shared_ptr<hierarchicalformfunction>(new h1pyramid(h1ddim));
        }
    }
    
    
    // If we arrive here it means the form function type name was incorrect:
    std::cout << "Error in 'selector' namespace: unknown form function type name '" << formfunctiontypename << "'" << std::endl;
    abort();
}
