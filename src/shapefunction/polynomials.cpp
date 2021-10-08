#include "polynomials.h"


polynomials::polynomials(std::vector<polynomial> input)
{
    mynumpolys = input.size();
    
    // Get the ki, eta and phi length of the polynomial encompassing them all:
    for (int pol = 0; pol < mynumpolys; pol++)
    {
        int curkilen = input[pol].mycoefficients.size();
        if (mykilen < curkilen)
            mykilen = curkilen;
        
        for (int k = 0; k < curkilen; k++)
        {
            int curetalen = input[pol].mycoefficients[k].size();
            if (myetalen < curetalen)
                myetalen = curetalen;
            
            for (int e = 0; e < curetalen; e++)
            {
                int curphilen = input[pol].mycoefficients[k][e].size();
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
        for (int k = 0; k < input[pol].mycoefficients.size(); k++)
        {
            for (int e = 0; e < input[pol].mycoefficients[k].size(); e++)
            {
                for (int p = 0; p < input[pol].mycoefficients[k][e].size(); p++)
                    mycoeffs[pol*mynummonomials+k*myetalen*myphilen+e*myphilen+p] = input[pol].mycoefficients[k][e][p];
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

void polynomials::evalatsingle(const std::vector<double>& evaluationpoint, int num, std::vector<double>& evaled)
{
    double ki = evaluationpoint[0], eta = evaluationpoint[1], phi = evaluationpoint[2];

    // Create the monomial value vector:
    std::vector<double> monomialval(mynummonomials);
    std::vector<int> kietaphipowers(3*mynummonomials);
    evaled = std::vector<double>((num+1)*mynumpolys,0);
    
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
                kietaphipowers[3*index+0] = k;
                kietaphipowers[3*index+1] = e;
                kietaphipowers[3*index+2] = p;
                index++;
                
                c *= phi;
            }
            b *= eta;
        }
        a *= ki;
    }

    // Multiply the coefficients by the monomial values:
    for (int pol = 0; pol < mynumpolys; pol++)
    {
        int p = pol*(num+1);
        for (int i = 0; i < mynummonomials; i++)
        {
            double cc = mycoeffs[pol*mynummonomials+i];
        
            evaled[p+0] += cc * monomialval[i];
            if (num > 0 && kietaphipowers[3*i+0] > 0)
                evaled[p+1] += cc * kietaphipowers[3*i+0] * monomialval[i-myetalen*myphilen]; // dki
            if (num > 1 && kietaphipowers[3*i+1] > 0)
                evaled[p+2] += cc * kietaphipowers[3*i+1] * monomialval[i-myphilen]; // deta
            if (num > 2 && kietaphipowers[3*i+2] > 0)
                evaled[p+3] += cc * kietaphipowers[3*i+2] * monomialval[i-1]; // dphi
        }
    }
}

polynomials polynomials::sum(std::vector<double>& weights)
{
    int num = weights.size()/mynumpolys;
    
    polynomials output;
    
    output.mynumpolys = num;
    
    output.mykilen = mykilen;
    output.myetalen = myetalen;
    output.myphilen = myphilen;
    
    output.mynummonomials = mynummonomials;
    
    output.mycoeffs = std::vector<double>(num * mynummonomials, 0);
    
    for (int n = 0; n < num; n++)
    {
        for (int p = 0; p < mynumpolys; p++)
        {
            double w = weights[n*mynumpolys+p];
            for (int c = 0; c < mynummonomials; c++)
                output.mycoeffs[n*mynummonomials+c] += w * mycoeffs[p*mynummonomials+c];
        }
    }
    
    return output;
}

void polynomials::print(void)
{
    std::cout << "Number of polynomials is " << mynumpolys << " with " << mynummonomials << " monomials (" << mykilen << "x" << myetalen << "x" << myphilen << ")" << std::endl;
    std::cout << "Coefficients are:" << std::endl;
    for (int i = 0; i < mycoeffs.size(); i++)
        std::cout << mycoeffs[i] << " ";
    std::cout << std::endl;
}

