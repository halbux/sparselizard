// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


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

#include "oneconstant.h"


namespace selector
{
    // Get a pointer to the selected form function:
    std::shared_ptr<hierarchicalformfunction> select(int elementtypenumber, std::string formfunctiontypename);
};

#endif
