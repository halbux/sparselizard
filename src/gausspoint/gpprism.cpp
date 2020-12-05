#include "gpprism.h"

int gpprism::count(int integrationorder)
{
    int numgplin = gpline::count(integrationorder);
    int numgptri = gptriangle::count(integrationorder);
    if (numgplin == -1 || numgptri == -1)
        return -1;
    else
        return numgplin*numgptri;
}

void gpprism::set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights)
{
    // General rule - based on the line and the triangle Gauss points:
    gausspoints gausspointsline(1, integrationorder);
    int numberoflinegausspoints = gausspointsline.count();
    std::vector<double> coordinatesline = gausspointsline.getcoordinates();
    std::vector<double> weightsline = gausspointsline.getweights();
    
    gausspoints gausspointstriangle(2, integrationorder);
    int numberoftrianglegausspoints = gausspointstriangle.count();
    std::vector<double> coordinatestriangle = gausspointstriangle.getcoordinates();
    std::vector<double> weightstriangle = gausspointstriangle.getweights();
    
    coordinates.resize(3*numberoflinegausspoints*numberoftrianglegausspoints);
    weights.resize(numberoflinegausspoints*numberoftrianglegausspoints);
    
    // Loop on all line Gauss points:
    int gp = 0;
    for (int i = 0; i < numberoflinegausspoints; i++)
    {
        // Loop on all triangle Gauss points:
        for (int j = 0; j < numberoftrianglegausspoints; j++)
        {
            coordinates[3*gp+0] = coordinatestriangle[3*j+0];
            coordinates[3*gp+1] = coordinatestriangle[3*j+1];
            coordinates[3*gp+2] = coordinatesline[3*i];
            weights[gp] = weightsline[i]*weightstriangle[j];
            
            gp++;
        }
    }
}




