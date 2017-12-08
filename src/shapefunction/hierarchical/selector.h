// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef SELECTOR_H
#define SELECTOR_H

#include <iostream>
#include <string>
#include <memory>

#include "h1point.h"
#include "h1line.h"
#include "h1triangle.h"
#include "h1quadrangle.h"
#include "h1tetrahedron.h"
#include "h1hexahedron.h"
#include "h1prism.h"
#include "h1pyramid.h"

#include "hcurlpoint.h"
#include "hcurlline.h"
#include "hcurltriangle.h"
#include "hcurlquadrangle.h"
#include "hcurltetrahedron.h"
#include "hcurlhexahedron.h"
#include "hcurlprism.h"
#include "hcurlpyramid.h"

#include "q6.h"
#include "h11.h"

#include "oneconstant.h"


namespace selector
{
    // Get a pointer to the selected form function:
    std::shared_ptr<hierarchicalformfunction> select(int elementtypenumber, std::string formfunctiontypename);
};

#endif
