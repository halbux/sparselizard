#include "lagrangepoint.h"


std::vector<double> lagrangepoint::getnodecoordinates([[maybe_unused]] int order)
{
    return std::vector<double> {0, 0, 0};
}

std::vector<polynomial> lagrangepoint::getformfunctionpolynomials(int order) // FIXME: use order
{
    std::vector<polynomial> formfunctionpoly(1);
    formfunctionpoly[0].set({{{1.0}}});

     return formfunctionpoly;
}
