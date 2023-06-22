#include "portrelation.h"

portrelation::portrelation(expression prtrel)
{
    prtrel.expand();

    // Extract the ports from the expression:
    std::vector<port> prts; std::vector<int> dtorders;
    prtrel.extractport(prts, dtorders, mycoefs, mynoportterm);

    int numports = prts.size();
    if (numports == 0)
    {
        logs log;
        log.msg() << "Error in 'portrelation' object: could not find any port in the port relation provided" << std::endl;
        log.error();
    }
    
    double pi = sl::getpi();
    
    // Get the harmonics in each port and the max harmonic found:
    int maxharm = 0;
    std::vector<std::vector<int>> portharms(numports);
    for (int p = 0; p < numports; p++)
    {
        portharms[p] = prts[p].getharmonics();
        int cmh = *max_element(portharms[p].begin(), portharms[p].end());
        maxharm = std::max(maxharm, cmh);
    }
    
    // Preallocate to one more for any possible dt(sin):
    myrawports.resize(maxharm+2);
    mycoefinds.resize(maxharm+2);
    mykcm.resize(maxharm+2);
    myfactors.resize(maxharm+2);

    // Extract the harmonics of each port and dispatch them to the appropriate relation:
    for (int p = 0; p < numports; p++)
    {
        int curdtorder = dtorders[p]; // max allowed is 2
        bool ismhport = prts[p].getpointer()->ismultiharmonic();
        
        for (int h = 0; h < portharms[p].size(); h++)
        {
            int harm = portharms[p][h];
            int harmfreq = harmonic::getfrequency(harm);
            std::shared_ptr<rawport> crp = prts[p].getpointer()->harmonic(harm);
            
            int targetharm = harm;
            int targetkcm = 0;
            double targetfactor = 1.0;
            int targetf0pow = 0;
            
            if (ismhport)
            {
                if (curdtorder == 1)
                {
                    if (harmonic::issine(harm))
                    {
                        targetharm = harm+1;
                        targetfactor = 2.0*pi*harmfreq;
                        targetf0pow = 1;
                    }
                    else
                    {
                        targetharm = harm-1;
                        targetfactor = -2.0*pi*harmfreq;
                        targetf0pow = 1;
                    }
                }
                if (curdtorder == 2)
                {
                    targetfactor = -4.0*pi*pi*harmfreq*harmfreq;
                    targetf0pow = 2;
                }
            }
            else
                targetkcm = curdtorder;

            if (targetharm > 0)
            {
                myrawports[targetharm].push_back(crp);
                mycoefinds[targetharm].push_back(p);
                mykcm[targetharm].push_back(targetkcm);
                myfactors[targetharm].push_back(std::make_pair(targetfactor, targetf0pow));
            }
        }
    }
}

int portrelation::count(void)
{
    int cnt = 0;
    for (int r = 1; r < myrawports.size(); r++)
    {
        if (r == 1 && mynoportterm.size() > 0)
        {
            cnt++;
            continue;
        }
    
        if (myrawports[r].size() > 0)
            cnt++;
    }
    
    return cnt;
}

std::vector<std::shared_ptr<rawport>> portrelation::getrawports(void)
{
    int numrp = 0;
    for (int r = 1; r < myrawports.size(); r++)
        numrp += myrawports[r].size();
        
    std::vector<std::shared_ptr<rawport>> output(numrp);
    
    int index = 0;
    for (int r = 1; r < myrawports.size(); r++)
    {
        for (int i = 0; i < myrawports[r].size(); i++)
            output[index+i] = myrawports[r][i];
        index += myrawports[r].size();
    }
    
    return output;
}

bool portrelation::hasnoportterm(void)
{
    return (mynoportterm.size() > 0);
}

double portrelation::evalnoportterm(void)
{
    return mynoportterm[0].evaluate();
}

void portrelation::evalrelations(int KCM, std::vector<std::shared_ptr<rawport>>& rps, std::vector<int>& relinds, std::vector<double>& relvals)
{
    // Count the number of terms:
    int numterms = 0;
    for (int r = 1; r < mykcm.size(); r++)
    {
        for (int i = 0; i < mykcm[r].size(); i++)
        {
            if (mykcm[r][i] == KCM)
                numterms++;
        }
    }
    
    // Preallocate:
    rps.resize(numterms);
    relinds.resize(numterms);
    relvals.resize(numterms);

    if (numterms == 0)
        return;
        
    // Calculate all coefficients:
    std::vector<double> coefvals(mycoefs.size());
    for (int i = 0; i < mycoefs.size(); i++)
        coefvals[i] = mycoefs[i].evaluate();

    // Populate:
    int ri = 0;
    int index = 0;
    
    for (int r = 1; r < mykcm.size(); r++)
    {
        for (int i = 0; i < mykcm[r].size(); i++)
        {
            if (mykcm[r][i] == KCM)
            {
                rps[index] = myrawports[r][i];
                relinds[index] = ri;
                
                // Calculate the coefficient:
                relvals[index] = coefvals[mycoefinds[r][i]] * myfactors[r][i].first;
                int f0pow = myfactors[r][i].second;
                if (f0pow != 0)
                    relvals[index] *= std::pow(universe::getfundamentalfrequency(), f0pow);
                
                index++;
            }
        }
        if (mykcm[r].size() > 0 || r == 1 && mynoportterm.size() > 0)
            ri++;
    }
}

