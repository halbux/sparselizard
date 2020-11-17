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
        if (val.size() < (size_t) h+1)
            val.resize(h+1);
        if (val[h].size() < 4)
            val[h].resize(4);
        if (val[h][i].size() < (size_t) j+1)
            val[h][i].resize(j+1);
        if (val[h][i][j].size() < (size_t) k+1)
            val[h][i][j].resize(k+1);
        if (val[h][i][j][k].size() < (size_t) l+1)
            val[h][i][j][k].resize(l+1);
        if (val[h][i][j][k][l].size() < 4)
            val[h][i][j][k][l].resize(4);
        if (val[h][i][j][k][l][m].size() < (size_t) n+1)
            val[h][i][j][k][l][m].resize(n+1);
    }
    
    if (ffpoly.size() < (size_t) h+1)
        ffpoly.resize(h+1);
    if (ffpoly[h].size() < 4)
        ffpoly[h].resize(4);
    if (ffpoly[h][i].size() < (size_t) j+1)
        ffpoly[h][i].resize(j+1);
    if (ffpoly[h][i][j].size() < (size_t) k+1)
        ffpoly[h][i][j].resize(k+1);
    if (ffpoly[h][i][j][k].size() < (size_t) l+1)
        ffpoly[h][i][j][k].resize(l+1);
    if (ffpoly[h][i][j][k][l].size() < (size_t) n+1)
        ffpoly[h][i][j][k][l].resize(n+1);

    // Add the polynomial:
    ffpoly[h][i][j][k][l][n] = poly;
}

void hierarchicalformfunctioncontainer::evaluate(std::vector<double> evaluationpoints)
{
    myevaluationpoints = evaluationpoints;

    for (size_t h = 0; h < val.size(); h++)
    {
        for (size_t i = 0; i < val[h].size(); i++)
        {
            for (size_t j = 0; j < val[h][i].size(); j++)
            {
                for (size_t k = 0; k < val[h][i][j].size(); k++)
                {
                    for (size_t l = 0; l < val[h][i][j][k].size(); l++)
                    {
                        for (size_t m = 0; m < val[h][i][j][k][l].size(); m++)
                        {
                            for (size_t n = 0; n < val[h][i][j][k][l][m].size(); n++)
                                val[h][i][j][k][l][m][n] = ffpoly[h][i][j][k][l][n].evalat(evaluationpoints, m);
                        }
                    }
                }
            }
        }
    }
}

densematrix hierarchicalformfunctioncontainer::tomatrix(int totalorientation, int order, int whichderivative, int component)
{
    std::vector<int> edgesorientations = orientation::getedgesorientationsfromtotalorientation(totalorientation, myelementtypenumber);
    std::vector<int> facesorientations = orientation::getfacesorientationsfromtotalorientation(totalorientation, myelementtypenumber);

    // Use an iterator to iterate through all form functions in the right order.
    hierarchicalformfunctioniterator myiterator(myformfunctiontypename, myelementtypenumber, order);

    // Create the numformfunc x numevalpoints dense matrix to output:
    densematrix valmat(myiterator.count(), myevaluationpoints.size()/3);
    
    for (int ff = 0; ff < myiterator.count(); ff++)
    {
        int h = myiterator.getformfunctionorder();
        int i = myiterator.getdimension();
        int j = myiterator.getnodeedgefacevolumeindex();
        int l = myiterator.getformfunctionindexincurrentorderinnodeedgefacevolume();

        // Get the orientation of the current node/edge/face/volume. 
        // Orientation is non zero only for edges and faces.
        int orientation = 0;
        if (i == 1)
            orientation = edgesorientations[j];
        if (i == 2)
            orientation = facesorientations[j];

        valmat.setrow(ff, val[h][i][j][orientation][l][whichderivative][component]);
        
        myiterator.next();
    }
    
    return valmat;
}

densematrix hierarchicalformfunctioncontainer::tomatrix(int h, int i, int j, int k, int l, int m, int n)
{
    return densematrix(1, myevaluationpoints.size()/3, val[h][i][j][k][l][m][n]);
}

void hierarchicalformfunctioncontainer::print(bool printallderivatives)
{
    std::cout.precision(17);
    for (size_t h = 0; h < val.size(); h++)
    {
        for (size_t i = 0; i < val[h].size(); i++)
        {
            for (size_t j = 0; j < val[h][i].size(); j++)
            {
                for (size_t k = 0; k < val[h][i][j].size(); k++)
                {
                    for (size_t l = 0; l < val[h][i][j][k].size(); l++)
                    {
                        for (size_t m = 0; m < val[h][i][j][k][l].size(); m++)
                        {
                            if (not(printallderivatives) && m > 0)
                                continue;
                            
                            for (size_t n = 0; n < val[h][i][j][k][l][m].size(); n++)
                            {
                                std::cout << "order " << h << "; dim " << i << "; sub elem index " << j << "; orientation " << k << "; number " << l << "; derivative " << m << "; comp " << n << ": ";
                                for (size_t o = 0; o < val[h][i][j][k][l][m][n].size(); o++)
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




