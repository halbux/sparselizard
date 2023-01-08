#include "hierarchicalformfunctioncontainer.h"
#include "hierarchicalformfunction.h"
#include "hierarchicalformfunctioniterator.h"

hierarchicalformfunctioncontainer::hierarchicalformfunctioncontainer(std::string formfunctiontypename, int elementtypenumber)
{
    myformfunctiontypename = formfunctiontypename;
    myelementtypenumber = elementtypenumber;
}

void hierarchicalformfunctioncontainer::set(int h, int i, int j, int k, int l, int n, polynomial& poly)
{   
    // Preallocate:
    for (int m = 0; m < 4; m++)
    {
        if (val.size() < h+1)
            val.resize(h+1);
        if (val[h].size() < 4)
            val[h].resize(4);
        if (val[h][i].size() < j+1)
            val[h][i].resize(j+1);
        if (val[h][i][j].size() < k+1)
            val[h][i][j].resize(k+1);
        if (val[h][i][j][k].size() < l+1)
            val[h][i][j][k].resize(l+1);
        if (val[h][i][j][k][l].size() < 4)
            val[h][i][j][k][l].resize(4);
        if (val[h][i][j][k][l][m].size() < n+1)
            val[h][i][j][k][l][m].resize(n+1);
    }
    
    if (ffpoly.size() < h+1)
        ffpoly.resize(h+1);
    if (ffpoly[h].size() < 4)
        ffpoly[h].resize(4);
    if (ffpoly[h][i].size() < j+1)
        ffpoly[h][i].resize(j+1);
    if (ffpoly[h][i][j].size() < k+1)
        ffpoly[h][i][j].resize(k+1);
    if (ffpoly[h][i][j][k].size() < l+1)
        ffpoly[h][i][j][k].resize(l+1);
    if (ffpoly[h][i][j][k][l].size() < n+1)
        ffpoly[h][i][j][k][l].resize(n+1);

    // Add the polynomial:
    ffpoly[h][i][j][k][l][n] = poly;
}

void hierarchicalformfunctioncontainer::evaluate(std::vector<double> evaluationpoints)
{
    if (myevaluationpoints == evaluationpoints)
        return;
        
    myevaluationpoints = evaluationpoints;

    for (int h = 0; h < val.size(); h++)
    {
        for (int i = 0; i < val[h].size(); i++)
        {
            for (int j = 0; j < val[h][i].size(); j++)
            {
                for (int k = 0; k < val[h][i][j].size(); k++)
                {
                    for (int l = 0; l < val[h][i][j][k].size(); l++)
                    {
                        for (int m = 0; m < val[h][i][j][k][l].size(); m++)
                        {
                            for (int n = 0; n < val[h][i][j][k][l][m].size(); n++)
                                val[h][i][j][k][l][m][n] = ffpoly[h][i][j][k][l][n].evalat(evaluationpoints, m);
                        }
                    }
                }
            }
        }
    }
}

densemat hierarchicalformfunctioncontainer::tomatrix(int totalorientation, int order, int whichderivative, int component)
{
    std::vector<int> edgesorientations = orientation::getedgesorientationsfromtotalorientation(totalorientation, myelementtypenumber);
    std::vector<int> facesorientations = orientation::getfacesorientationsfromtotalorientation(totalorientation, myelementtypenumber);

    // Use an iterator to iterate through all form functions in the right order.
    hierarchicalformfunctioniterator myiterator(myformfunctiontypename, myelementtypenumber, order);

    // Create the numformfunc x numevalpoints dense matrix to output:
    densemat valmat(myiterator.count(), myevaluationpoints.size()/3);
    
    for (int ff = 0; ff < myiterator.count(); ff++)
    {
        int h = myiterator.getformfunctionorder();
        int i = myiterator.getdimension();
        int j = myiterator.getnodeedgefacevolumeindex();
        int l = myiterator.getformfunctionindexincurrentorderinnodeedgefacevolume();

        // Get the orientation of the current node/edge/face/volume. 
        // Orientation is non zero only for edges and faces.
        int orientation = 0;
        if (val[h][i][j].size() > 1) // e.g. not for h1d
        {
            if (i == 1)
                orientation = edgesorientations[j];
            if (i == 2)
                orientation = facesorientations[j];
        }

        valmat.setrow(ff, val[h][i][j][orientation][l][whichderivative][component]);
        
        myiterator.next();
    }
    
    return valmat;
}

densemat hierarchicalformfunctioncontainer::tomatrix(int h, int i, int j, int k, int l, int m, int n)
{
    return densemat(1, myevaluationpoints.size()/3, val[h][i][j][k][l][m][n]);
}

void hierarchicalformfunctioncontainer::print(bool printallderivatives)
{
    std::cout.precision(17);
    for (int h = 0; h < val.size(); h++)
    {
        for (int i = 0; i < val[h].size(); i++)
        {
            for (int j = 0; j < val[h][i].size(); j++)
            {
                for (int k = 0; k < val[h][i][j].size(); k++)
                {
                    for (int l = 0; l < val[h][i][j][k].size(); l++)
                    {
                        for (int m = 0; m < val[h][i][j][k][l].size(); m++)
                        {
                            if (not(printallderivatives) && m > 0)
                                continue;
                            
                            for (int n = 0; n < val[h][i][j][k][l][m].size(); n++)
                            {
                                std::cout << "order " << h << "; dim " << i << "; sub elem index " << j << "; orientation " << k << "; number " << l << "; derivative " << m << "; comp " << n << ": ";
                                for (int o = 0; o < val[h][i][j][k][l][m][n].size(); o++)
                                    std::cout << val[h][i][j][k][l][m][n][o] << " ";
                                std::cout << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

