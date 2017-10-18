#include "lagrangepyramid.h"


std::vector<double> lagrangepyramid::getnodecoordinates(int order)
{
     switch (order)
     {
          case 1:
               return std::vector<double> {-1.0, -1.0, 0, 1.0, -1.0, 0, 1.0, 1.0, 0, -1.0, 1.0, 0.0, 0.0, 0.0, 1.0};
         default:
               std::cout << "Error in 'lagrangepyramid' namespace: coordinates of order 2 and above not defined" << std::endl;
               abort();
               break;
     }
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

                    // In Zaglmayr's thesis the reference elements are shifted and deformed.
                    // Variable change to correspond to our reference element definition:
                    // ki := (ki+1)/2; eta := (eta+1)/2;
                    ki = 0.5*(ki+1); eta = 0.5*(eta+1);
              
                    formfunctionpoly[0] = (1.0-ki)*(1.0-eta)*(1.0-phi);
                    formfunctionpoly[1] = ki*(1.0-eta)*(1.0-phi);
                    formfunctionpoly[2] = ki*eta*(1.0-phi);
                    formfunctionpoly[3] = (1.0-ki)*eta*(1.0-phi);
                    formfunctionpoly[4] = phi;
              
//                formfunctionpoly[0].set({{{0.5, -0.5}, {-0.5, 0.5}}, {{-0.5, 0.5}}});
//                formfunctionpoly[1].set({{{0.0, 0.0}, {0.0, 0.0}}, {{0.5, -0.5}}});
//                formfunctionpoly[2].set({{{0.0, 0.0}, {0.5, -0.5}}, {{0.0, 0.0}}});
//                formfunctionpoly[3].set({{{0.5, 0.5}, {-0.5, -0.5}}, {{-0.5, -0.5}}});
//                formfunctionpoly[4].set({{{0.0, 0.0}, {0.0, 0.0}}, {{0.5, 0.5}}});
               break;
          }
          default:
               std::cout << "Error in 'lagrangepyramid' namespace: Lagrange form functions of order 2 and above not defined" << std::endl;
               abort();
               break;
     }

     return formfunctionpoly;
}