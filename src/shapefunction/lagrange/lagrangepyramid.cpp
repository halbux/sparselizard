#include "lagrangepyramid.h"


std::vector<double> lagrangepyramid::getnodecoordinates(int order)
{
     switch (order)
     {
          case 1:
               return std::vector<double> {-1.0, -1.0, 0, 1.0, -1.0, 0, 1.0, 1.0, 0, -1.0, 1.0, 0.0, 0.0, 0.0, 1.0};
         default:
               logs log;
               log.msg() << "Error in 'lagrangepyramid' namespace: coordinates of order 2 and above not defined" << std::endl;
               log.error();
               break;
     }
     
     throw std::runtime_error(""); // fix return warning
}

std::vector<polynomial> lagrangepyramid::getformfunctionpolynomials(int order)
{
     element pyramid(7,order);
     std::vector<polynomial> formfunctionpoly(pyramid.countcurvednodes());

     switch (order)
     {
          case 1:
          {
                polynomial ki, eta, phi;
                ki.set({{{}},{{{1.0}}}});
                eta.set({{{},{1.0}}});
                phi.set({{{0.0,1.0}}});

                formfunctionpoly[0] = 0.25*(1.0-ki)*(1.0-eta)*(1.0-phi);
                formfunctionpoly[1] = 0.25*(1.0+ki)*(1.0-eta)*(1.0-phi);
                formfunctionpoly[2] = 0.25*(1.0+ki)*(1.0+eta)*(1.0-phi);
                formfunctionpoly[3] = 0.25*(1.0-ki)*(1.0+eta)*(1.0-phi);
                formfunctionpoly[4] = phi;
                break;
          }
          default:
               logs log;
               log.msg() << "Error in 'lagrangepyramid' namespace: Lagrange form functions of order 2 and above not defined" << std::endl;
               log.error();
               break;
     }

     return formfunctionpoly;
}
