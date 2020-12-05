#include "gphexahedron.h"

int gphexahedron::count(int integrationorder)
{
    int numgplin = gpline::count(integrationorder);
    if (numgplin == -1)
        return -1;
    else
        return numgplin*numgplin*numgplin;
}

void gphexahedron::set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights)
{
    // General rule - based on the 1D Gauss points:
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
                coordinates[3*gp+0] = coordinatesline[3*i];
                coordinates[3*gp+1] = coordinatesline[3*j];
                coordinates[3*gp+2] = coordinatesline[3*k];
                weights[gp] = weightsline[i]*weightsline[j]*weightsline[k];

                gp++;
            }
        }
    }
}




