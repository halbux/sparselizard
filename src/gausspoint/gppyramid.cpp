#include "gppyramid.h"

void gppyramid::set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights)
{
//     std::cout << "Error 'in gppyramid' namespace: Gauss points have not been defined yet for pyramids" << std::endl;
//     abort();
    
    // General rule - based on the Duffy transformation of a hexahedron:
    gausspoints gausspointsline(1, integrationorder);
    int numberoflinegausspoints = gausspointsline.count();
    std::vector<double> coordinatesline = gausspointsline.getcoordinates();
    std::vector<double> weightsline = gausspointsline.getweights();

    coordinates.resize(3*pow(numberoflinegausspoints,3));
    weights.resize(pow(numberoflinegausspoints,3));
    
    int gp = 0;
    for (int i = 0; i < numberoflinegausspoints; i++)
    {
        for (int j = 0; j < numberoflinegausspoints; j++)
        {
            for (int k = 0; k < numberoflinegausspoints; k++)
            {
                // For the pyramid the phi coordinates ranges from 0 to 1 while for 
                // the hexahedron it ranges from -1 to 1 --> scaling factor:
                coordinates[3*gp+2] = 0.5*(coordinatesline[3*k] + 1);
                
                coordinates[3*gp+0] = coordinatesline[3*i]*(1-coordinates[3*gp+2]);
                coordinates[3*gp+1] = coordinatesline[3*j]*(1-coordinates[3*gp+2]);

                weights[gp] = weightsline[i]*weightsline[j]*weightsline[k];
                // Volume of hex is 8, of pyramid is 4/3.
                weights[gp] *= 1.0/8.0 * pow(4.0/3.0,3);
                
                gp++; //////// WRONG.... :(
            }
        }
    }
    
}




