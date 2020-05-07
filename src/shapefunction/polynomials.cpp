#include "polynomials.h"


polynomials::polynomials(std::vector<polynomial> input)
{
    mypolys = input;
    mynumpolys = mypolys.size();
    
    // Get the ki, eta and phi length of the polynomial encompassing them all:
    for (int pol = 0; pol < mynumpolys; pol++)
    {
        int curkilen = mypolys[pol].mycoefficients.size();
        if (mykilen < curkilen)
            mykilen = curkilen;
        
        for (int k = 0; k < curkilen; k++)
        {
            int curetalen = mypolys[pol].mycoefficients[k].size();
            if (myetalen < curetalen)
                myetalen = curetalen;
            
            for (int e = 0; e < curetalen; e++)
            {
                int curphilen = mypolys[pol].mycoefficients[k][e].size();
                if (myphilen < curphilen)
                    myphilen = curphilen;
            }
        }
    }
    
    // Populate the coefficients container:
    mynummonomials = mykilen*myetalen*myphilen;
    mycoeffs = std::vector<double>(mynumpolys * mynummonomials, 0.0);
    for (int pol = 0; pol < mynumpolys; pol++)
    {
        for (int k = 0; k < mypolys[pol].mycoefficients.size(); k++)
        {
            for (int e = 0; e < mypolys[pol].mycoefficients[k].size(); e++)
            {
                for (int p = 0; p < mypolys[pol].mycoefficients[k][e].size(); p++)
                    mycoeffs[pol*mynummonomials+k*myetalen*myphilen+e*myphilen+p] = mypolys[pol].mycoefficients[k][e][p];
            }
        }
    }
}

void polynomials::evalatsingle(const std::vector<double>& evaluationpoint, std::vector<double>& evaled)
{
    double ki = evaluationpoint[0], eta = evaluationpoint[1], phi = evaluationpoint[2];


    // Create the monomial value vector:
    std::vector<double> monomialval(mynummonomials);
    evaled = std::vector<double>(mynumpolys,0);
    
    int index = 0;
    double a = 1, b = 1, c = 1;
    for (int k = 0; k < mykilen; k++)
    {   
        b = a;
        for (int e = 0; e < myetalen; e++)
        {    
            c = b;
            for (int p = 0; p < myphilen; p++)
            {
                monomialval[index] = c;   
                index++;
                
                c *= phi;
            }
            b *= eta;
        }
        a *= ki;
    }

    // Multiply the coefficients by the monomial values:
    index = 0;
    for (int pol = 0; pol < mynumpolys; pol++)
    {
        for (int i = 0; i < mynummonomials; i++)
        {
            evaled[pol] += mycoeffs[index]*monomialval[i];
            index++;
        }
    }
}

void polynomials::print(void)
{
    std::cout << "Number of polynomials is " << mynumpolys << " with " << mynummonomials << " monomials (" << mykilen << "x" << myetalen << "x" << myphilen << ")" << std::endl;
    std::cout << "Coefficients are:" << std::endl;
    for (int i = 0; i < mycoeffs.size(); i++)
        std::cout << mycoeffs[i] << " ";
    std::cout << std::endl;
}

