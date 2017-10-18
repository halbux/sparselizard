#include "selector.h"


std::shared_ptr<hierarchicalformfunction> selector::select(int elementtypenumber, std::string formfunctiontypename)
{
    if (formfunctiontypename.compare("h1") == 0)
	{
        switch (elementtypenumber)
        {
            case 0:
                return std::shared_ptr<hierarchicalformfunction>(new h1point);
            case 1:
                return std::shared_ptr<hierarchicalformfunction>(new h1line);
            case 2:
                return std::shared_ptr<hierarchicalformfunction>(new h1triangle);
            case 3:
                return std::shared_ptr<hierarchicalformfunction>(new h1quadrangle);
            case 4:
                return std::shared_ptr<hierarchicalformfunction>(new h1tetrahedron);
            case 5:
                return std::shared_ptr<hierarchicalformfunction>(new h1hexahedron);
            case 6:
                return std::shared_ptr<hierarchicalformfunction>(new h1prism);
            case 7:
                return std::shared_ptr<hierarchicalformfunction>(new h1pyramid);
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
    
    if (formfunctiontypename.compare("q6") == 0)
	{
        switch (elementtypenumber)
        {
            case 0:
                return std::shared_ptr<hierarchicalformfunction>(new h1point);
            case 1:
                return std::shared_ptr<hierarchicalformfunction>(new h1line);
            case 2:
                return std::shared_ptr<hierarchicalformfunction>(new h1triangle);
            case 3:
                return std::shared_ptr<hierarchicalformfunction>(new q6);
            case 4:
                return std::shared_ptr<hierarchicalformfunction>(new h1tetrahedron);
            case 5:
                return std::shared_ptr<hierarchicalformfunction>(new h1hexahedron);
            case 6:
                return std::shared_ptr<hierarchicalformfunction>(new h1prism);
            case 7:
                return std::shared_ptr<hierarchicalformfunction>(new h1pyramid);
        }
	}
    
    if (formfunctiontypename.compare("h11") == 0)
	{
        switch (elementtypenumber)
        {
            case 0:
                return std::shared_ptr<hierarchicalformfunction>(new h1point);
            case 1:
                return std::shared_ptr<hierarchicalformfunction>(new h1line);
            case 2:
                return std::shared_ptr<hierarchicalformfunction>(new h1triangle);
            case 3:
                return std::shared_ptr<hierarchicalformfunction>(new h1quadrangle);
            case 4:
                return std::shared_ptr<hierarchicalformfunction>(new h1tetrahedron);
            case 5:
                return std::shared_ptr<hierarchicalformfunction>(new h11);
            case 6:
                return std::shared_ptr<hierarchicalformfunction>(new h1prism);
            case 7:
                return std::shared_ptr<hierarchicalformfunction>(new h1pyramid);
        }
	}
    
    
    // If we arrive here it means the form function type name was incorrect:
    std::cout << "Error in 'selector' namespace: unknown form function type name '" << formfunctiontypename << "'" << std::endl;
    abort();
}
