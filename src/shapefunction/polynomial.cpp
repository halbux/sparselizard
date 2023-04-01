#include "polynomial.h"


void polynomial::set(const std::vector<std::vector<std::vector<double>>>& coefficients)
{
    mycoefficients = coefficients;
}

void polynomial::print(void)
{
    for (int kiorder = 0; kiorder < mycoefficients.size(); kiorder++)
    {
        for (int etaorder = 0; etaorder < mycoefficients[kiorder].size(); etaorder++)
        {
            for (int phiorder = 0; phiorder < mycoefficients[kiorder][etaorder].size(); phiorder++)
            {
                double currentcoef = mycoefficients[kiorder][etaorder][phiorder];
                // Skip 0 coefficients:
                if (currentcoef != 0)
                {
                    // Print the coefficient:
                    if (currentcoef > 0)
                        std::cout << '+' << currentcoef;
                    else
                        std::cout << currentcoef;
                
                    // Print the ki, eta, phi product:
                    if (kiorder != 0)
                    {
                        if (kiorder == 1)
                            std::cout << "*ki";
                        else
                            std::cout << "*ki^" << kiorder;
                    }
                    if (etaorder != 0)
                    {
                        if (etaorder == 1)
                            std::cout << "*eta";
                        else
                            std::cout << "*eta^" << etaorder;
                    }
                    if (phiorder != 0)
                    {
                        if (phiorder == 1)
                            std::cout << "*phi";
                        else
                            std::cout << "*phi^" << phiorder;
                    }
                }
            }
        }
    }
    std::cout << std::endl;
}

std::vector<double> polynomial::evalat(const std::vector<double>& evaluationpoints, int whichderivative)
{
    int firstkiorder = 0;
    int firstetaorder = 0;
    int firstphiorder = 0;
    if (whichderivative == 1)
        firstkiorder = 1;
    if (whichderivative == 2)
        firstetaorder = 1;
    if (whichderivative == 3)
        firstphiorder = 1;

    int numberofevaluationpoints = evaluationpoints.size()/3;
    std::vector<double> val(numberofevaluationpoints,0);
    
    for (int i = 0; i < numberofevaluationpoints; i++)
    {
        // Current ki, eta and phi values:
        double ki = evaluationpoints[3*i+0];
        double eta = evaluationpoints[3*i+1];
        double phi = evaluationpoints[3*i+2];
    
        double kipower = 1;
    
        // Loop on all coefficients:
        for (int kiorder = firstkiorder; kiorder < mycoefficients.size(); kiorder++)
        {        
            double etapower = 1;
            for (int etaorder = firstetaorder; etaorder < mycoefficients[kiorder].size(); etaorder++)
            {
                double phipower = 1;
                for (int phiorder = firstphiorder; phiorder < mycoefficients[kiorder][etaorder].size(); phiorder++)
                {
                    // Skip zero coefficients:
                    if (mycoefficients[kiorder][etaorder][phiorder] != 0)
                    {
                        // In case a derivative has to be computed one has to take into account the extra factor:
                        if (whichderivative == 0)
                            val[i] += mycoefficients[kiorder][etaorder][phiorder] * kipower * etapower * phipower;
                        if (whichderivative == 1)
                            val[i] += mycoefficients[kiorder][etaorder][phiorder] * kipower * etapower * phipower * kiorder;
                        if (whichderivative == 2)
                            val[i] += mycoefficients[kiorder][etaorder][phiorder] * kipower * etapower * phipower * etaorder;
                        if (whichderivative == 3)
                            val[i] += mycoefficients[kiorder][etaorder][phiorder] * kipower * etapower * phipower * phiorder;
                    }
                    phipower = phi*phipower;
                }
                etapower = eta*etapower;
            }
            kipower = ki*kipower;
        }
    }
    return val;
}

polynomial polynomial::operator*(polynomial tomultiply)
{
    polynomial output({});
    
    // Loop on all coefficients of 'tomultiply' (descendingly to reduce the number of .resize calls on 'output):
    for (int kiorder = tomultiply.mycoefficients.size()-1; kiorder >= 0; kiorder--)
    {    
        for (int etaorder = tomultiply.mycoefficients[kiorder].size()-1; etaorder >= 0; etaorder--)
        {
            for (int phiorder = tomultiply.mycoefficients[kiorder][etaorder].size()-1; phiorder >= 0; phiorder--)
            {
                // If the current coefficient of 'tomultiply' is 0 skip it:
                if (tomultiply.mycoefficients[kiorder][etaorder][phiorder] == 0)
                    continue;
                
                // Loop on all coefficients of 'this':
                for (int kiorderthis = mycoefficients.size()-1; kiorderthis >= 0; kiorderthis--)
                {    
                    for (int etaorderthis = mycoefficients[kiorderthis].size()-1; etaorderthis >= 0; etaorderthis--)
                    {
                        for (int phiorderthis = mycoefficients[kiorderthis][etaorderthis].size()-1; phiorderthis >= 0; phiorderthis--)
                        {
                            // If the current coefficient of 'this' is 0 skip it:
                            if (mycoefficients[kiorderthis][etaorderthis][phiorderthis] == 0)
                                continue;
                                
                            // First preallocate 'output' for the current product term (if required):
                            if (output.mycoefficients.size() < kiorder+kiorderthis+1)
                                output.mycoefficients.resize(kiorder+kiorderthis+1);
                            if (output.mycoefficients[kiorder+kiorderthis].size() < etaorder+etaorderthis+1)
                                output.mycoefficients[kiorder+kiorderthis].resize(etaorder+etaorderthis+1);
                            if (output.mycoefficients[kiorder+kiorderthis][etaorder+etaorderthis].size() < phiorder+phiorderthis+1)
                                output.mycoefficients[kiorder+kiorderthis][etaorder+etaorderthis].resize(phiorder+phiorderthis+1);
                            
                            // Multiply and add to 'output':
                            output.mycoefficients[kiorder+kiorderthis][etaorder+etaorderthis][phiorder+phiorderthis] += mycoefficients[kiorderthis][etaorderthis][phiorderthis] * tomultiply.mycoefficients[kiorder][etaorder][phiorder];
                        }
                    }
                }
            }
        }
    }
    return output;
}

polynomial polynomial::operator+(polynomial toadd)
{
    polynomial output = *this;
    
    if (output.mycoefficients.size() < toadd.mycoefficients.size())
        output.mycoefficients.resize(toadd.mycoefficients.size());
            
    // Loop on all coefficients of 'toadd':
    for (int kiorder = 0; kiorder < toadd.mycoefficients.size(); kiorder++)
    {        
        for (int etaorder = 0; etaorder < toadd.mycoefficients[kiorder].size(); etaorder++)
        {
            if (output.mycoefficients[kiorder].size() < toadd.mycoefficients[kiorder].size())
                output.mycoefficients[kiorder].resize(toadd.mycoefficients[kiorder].size());
                
            for (int phiorder = 0; phiorder < toadd.mycoefficients[kiorder][etaorder].size(); phiorder++)
            {
                if (output.mycoefficients[kiorder][etaorder].size() < toadd.mycoefficients[kiorder][etaorder].size())
                    output.mycoefficients[kiorder][etaorder].resize(toadd.mycoefficients[kiorder][etaorder].size());
                // Add the two coefficients:
                output.mycoefficients[kiorder][etaorder][phiorder] += toadd.mycoefficients[kiorder][etaorder][phiorder];
            }
        }
    }
    return output;
}

polynomial polynomial::operator-(polynomial tosubtract)
{
    polynomial output = *this;
    
    if (output.mycoefficients.size() < tosubtract.mycoefficients.size())
        output.mycoefficients.resize(tosubtract.mycoefficients.size());
            
    // Loop on all coefficients of 'tosubtract':
    for (int kiorder = 0; kiorder < tosubtract.mycoefficients.size(); kiorder++)
    {        
        for (int etaorder = 0; etaorder < tosubtract.mycoefficients[kiorder].size(); etaorder++)
        {
            if (output.mycoefficients[kiorder].size() < tosubtract.mycoefficients[kiorder].size())
                output.mycoefficients[kiorder].resize(tosubtract.mycoefficients[kiorder].size());
                
            for (int phiorder = 0; phiorder < tosubtract.mycoefficients[kiorder][etaorder].size(); phiorder++)
            {
                if (output.mycoefficients[kiorder][etaorder].size() < tosubtract.mycoefficients[kiorder][etaorder].size())
                    output.mycoefficients[kiorder][etaorder].resize(tosubtract.mycoefficients[kiorder][etaorder].size());
                // Subtract the two coefficients:
                output.mycoefficients[kiorder][etaorder][phiorder] -= tosubtract.mycoefficients[kiorder][etaorder][phiorder];
            }
        }
    }
    return output;
}

polynomial polynomial::operator+()
{
    return *this;
}

polynomial polynomial::operator-()
{
    polynomial output = *this;
    
    // Loop on all coefficients:
    for (int kiorder = 0; kiorder < mycoefficients.size(); kiorder++)
    {
        for (int etaorder = 0; etaorder < mycoefficients[kiorder].size(); etaorder++)
        {
            for (int phiorder = 0; phiorder < mycoefficients[kiorder][etaorder].size(); phiorder++)
            {
                // Flip the sign of the coefficient:
                output.mycoefficients[kiorder][etaorder][phiorder] = -mycoefficients[kiorder][etaorder][phiorder];
            }
        }
    }
    return output;
}

polynomial polynomial::operator*(double val)
{
    polynomial output = *this;
    
    // Loop on all coefficients:
    for (int kiorder = 0; kiorder < mycoefficients.size(); kiorder++)
    {
        for (int etaorder = 0; etaorder < mycoefficients[kiorder].size(); etaorder++)
        {
            for (int phiorder = 0; phiorder < mycoefficients[kiorder][etaorder].size(); phiorder++)
            {
                // Multiply the coefficient by 'val':
                output.mycoefficients[kiorder][etaorder][phiorder] = val*mycoefficients[kiorder][etaorder][phiorder];
            }
        }
    }
    return output;
}

polynomial polynomial::operator+(double val)
{
    polynomial output = *this;
    // Make sure the [0][0][0] coefficient is defined:
    if (output.mycoefficients.size() < 1)
        output.mycoefficients.resize(1);
    if (output.mycoefficients[0].size() < 1)
        output.mycoefficients[0].resize(1);
    if (output.mycoefficients[0][0].size() < 1)
        output.mycoefficients[0][0].resize(1);
        
    output.mycoefficients[0][0][0] = output.mycoefficients[0][0][0] + val;
    return output;
}

polynomial polynomial::operator-(double val)
{
    polynomial output = *this;
    // Make sure the [0][0][0] coefficient is defined:
    if (output.mycoefficients.size() < 1)
        output.mycoefficients.resize(1);
    if (output.mycoefficients[0].size() < 1)
        output.mycoefficients[0].resize(1);
    if (output.mycoefficients[0][0].size() < 1)
        output.mycoefficients[0][0].resize(1);
        
    output.mycoefficients[0][0][0] = output.mycoefficients[0][0][0] - val;
    return output;
}

void polynomial::dki(void)
{
    // Shift all ki vectors back down one order - dump order 0 in ki:
    for (int kiorder = 1; kiorder < mycoefficients.size(); kiorder++)
        mycoefficients[kiorder-1] = mycoefficients[kiorder];
    if (mycoefficients.size() > 0)
        mycoefficients.pop_back();
    // One also has to take into account the extra factor coming from a derivative:    
    for (int kiorder = 1; kiorder < mycoefficients.size(); kiorder++)
    {
        for (int etaorder = 0; etaorder < mycoefficients[kiorder].size(); etaorder++)
        {
            for (int phiorder = 0; phiorder < mycoefficients[kiorder][etaorder].size(); phiorder++)
                mycoefficients[kiorder][etaorder][phiorder] = (kiorder+1) * mycoefficients[kiorder][etaorder][phiorder];
        }
    }
}

void polynomial::deta(void)
{
    // Shift all eta vectors back down one order - dump order 0 in eta:
    for (int kiorder = 0; kiorder < mycoefficients.size(); kiorder++)
    {
        for (int etaorder = 1; etaorder < mycoefficients[kiorder].size(); etaorder++)
            mycoefficients[kiorder][etaorder-1] = mycoefficients[kiorder][etaorder];
        if (mycoefficients[kiorder].size() > 0)
            mycoefficients[kiorder].pop_back();
        // One also has to take into account the extra factor coming from a derivative:    
        for (int etaorder = 1; etaorder < mycoefficients[kiorder].size(); etaorder++)
        {
            for (int phiorder = 0; phiorder < mycoefficients[kiorder][etaorder].size(); phiorder++)
                mycoefficients[kiorder][etaorder][phiorder] = (etaorder+1) * mycoefficients[kiorder][etaorder][phiorder];
        }
    }
}

void polynomial::dphi(void)
{
    // Shift all phi vectors back down one order - dump order 0 in phi:
    for (int kiorder = 0; kiorder < mycoefficients.size(); kiorder++)
    {
        for (int etaorder = 0; etaorder < mycoefficients[kiorder].size(); etaorder++)
        {
            for (int phiorder = 1; phiorder < mycoefficients[kiorder][etaorder].size(); phiorder++)
                mycoefficients[kiorder][etaorder][phiorder-1] = mycoefficients[kiorder][etaorder][phiorder];
            if (mycoefficients[kiorder][etaorder].size() > 0)
                mycoefficients[kiorder][etaorder].pop_back();
            // One also has to take into account the extra factor coming from a derivative:    
            for (int phiorder = 1; phiorder < mycoefficients[kiorder][etaorder].size(); phiorder++)
                mycoefficients[kiorder][etaorder][phiorder] = (phiorder+1) * mycoefficients[kiorder][etaorder][phiorder];
        }
    }
}

polynomial polynomial::derivative(int whichderivative)
{
    switch (whichderivative)
    {
        case 0:
        {
            polynomial dkipoly = *this;
            dkipoly.dki();
            return dkipoly;
        }
        case 1:
        {
            polynomial detapoly = *this;
            detapoly.deta();
            return detapoly;
        }
        case 2:
        {
            polynomial dphipoly = *this;
            dphipoly.dphi();
            return dphipoly;
        }
    }
    
    throw std::runtime_error(""); // fix return warning
}



polynomial operator*(double val, polynomial poly) {return poly*val;}
polynomial operator+(double val, polynomial poly) {return poly+val;}
polynomial operator-(double val, polynomial poly) {return -poly+val;}





