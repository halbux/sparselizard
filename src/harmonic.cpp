#include "harmonic.h"


int harmonic::getfrequency(int harmonicnumber)
{
    return (harmonicnumber - harmonicnumber%2 )/2;
}

bool harmonic::issine(int harmonicnumber)
{
    return (harmonicnumber%2 == 0);
}

bool harmonic::iscosine(int harmonicnumber)
{
    return (harmonicnumber%2 != 0);
}

int harmonic::getharmonicnumber(int frequency, bool issine)
{
    if (issine)
        return frequency*2;
    else
        return frequency*2+1;
}

int harmonic::getharmonicnumber(std::string input)
{
    if (input.size() >= 4 && (input.compare(0,3,"sin") == 0 || input.compare(0,3,"cos") == 0) && std::all_of(input.begin()+3, input.end(), ::isdigit) == true)
    {
        bool issine = true;
        if (input.compare(0,3,"cos") == 0)
            issine = false;
        
        input.erase(0,3);
        int freqindex = std::stoi(input);
        
        return getharmonicnumber(freqindex, issine);
    }
    logs log;
    log.msg() << "Error in 'harmonic' namespace: '" << input << "' is not a valid harmonic name (correct form is cos0, sin1, cos1, sin2, cos2, ...)" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<std::pair<int,double>> harmonic::getproduct(int harm1, int harm2)
{
    // Special case if any of the two inputs is cos0:
    if (harm1 == 1)
        return {std::make_pair(harm2,1)};
    if (harm2 == 1)
        return {std::make_pair(harm1,1)};
    
    // The cos(a) x sin(b) case is treated in sin(b) x cos(a) case:
    if (iscosine(harm1) && issine(harm2))
        std::swap(harm1,harm2);
        
    // If none is the constant harmonic:
    int freq1 = getfrequency(harm1);
    int freq2 = getfrequency(harm2);
    
    // sin(a) x sin(b) case:
    if (issine(harm1) && issine(harm2))
        return {std::make_pair(getharmonicnumber(std::abs(freq1-freq2), false),0.5), std::make_pair(getharmonicnumber(freq1+freq2, false),-0.5)};
    // cos(a) x cos(b) case:
    if (iscosine(harm1) && iscosine(harm2))
        return {std::make_pair(getharmonicnumber(std::abs(freq1-freq2), false),0.5), std::make_pair(getharmonicnumber(freq1+freq2, false),0.5)};
    // sin(a) x cos(b) case:
    if (issine(harm1) && iscosine(harm2))
    {
        if (freq1 != freq2)
        {
            if (freq1 > freq2)
                return {std::make_pair(getharmonicnumber(freq1+freq2, true),0.5), std::make_pair(getharmonicnumber(freq1-freq2, true),0.5)};
            else
                return {std::make_pair(getharmonicnumber(freq1+freq2, true),0.5), std::make_pair(getharmonicnumber(freq2-freq1, true),-0.5)};
        }
        else
            return {std::make_pair(getharmonicnumber(std::abs(freq1+freq2), true),0.5)};
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<std::pair<int,double>> harmonic::getproduct(int harm1, int harm2, int harm2timederivativeorder)
{
    if (harm2timederivativeorder == 0)
        return getproduct(harm1, harm2);
    
    // Derivating harmonic 1 gives zero:
    if (harm2 == 1)
        return {};

    // Call the classical 'getproduct' function but with modified harm 2.
    // For odd derivation orders if harm2 is on a sine it is put on a cosine and vice versa.
    int newharm2 = harm2;
    if (harm2timederivativeorder%2 == 1)
        newharm2 = newharm2 + 1 - 2*(newharm2%2);

    std::vector<std::pair<int,double>> harmsofproduct = getproduct(harm1, newharm2);

    // Multiply the coefficient by the extra 2*pi*fi like factor from the time derivation:
    for (int p = 0; p < harmsofproduct.size(); p++)
        harmsofproduct[p].second *= getderivationfactor(harm2timederivativeorder, harm2);

    return harmsofproduct;
}

double harmonic::getderivationfactor(int timederivativeorder, int harm)
{
    if (timederivativeorder > 0)
    {
        if (harm < 2)
            return 0;
    
        double pi = 3.141592653589793;
        double f0 = universe::getfundamentalfrequency();

        double mysign = 1;
        if (issine(harm))
        {
            for (int i = 2; i <= timederivativeorder; i = i+2)
                mysign *= -1;
        }
        else
        {
            for (int i = 1; i <= timederivativeorder; i = i+2)
                mysign *= -1;
        }
        return mysign * std::pow(2*pi*f0 * getfrequency(harm), timederivativeorder);
    }
    else
        return 1;
}

std::vector<std::vector<densemat>> harmonic::timederivative(int timederivativeorder, std::vector<std::vector<densemat>> input)
{
    if (timederivativeorder == 0)
        return input;
   
    for (int h = 0; h < input.size(); h++)
    {
        if (input[h].size() == 1)
            input[h][0].multiplyelementwise(getderivationfactor(timederivativeorder, h));
    }
    
    // Flip the sin and cos harmonics if required:
    if (timederivativeorder%2 == 1)
    {
        if (input.size()%2 != 0)
            input.push_back({});
        
        // Harmonic 1 is zero anyway. Start at 2.
        for (int h = 2; h < input.size(); h = h+2)
        {
            std::vector<densemat> temp = input[h];
            input[h] = input[h+1];
            input[h+1] = temp;
        }
    }
    return input;
}
