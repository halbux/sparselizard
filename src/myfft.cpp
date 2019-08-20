#include "myfft.h"


std::vector<std::vector<densematrix>> myfft::fft(densematrix input, long int mym, long int myn)
{
    // Number of time evaluations.
    long int numtimeevals = input.countrows();
    // Number of 1D transforms to perform:
    long int numtransforms = input.countcolumns();


    double pi = 3.141592653589793238;

    double* inputvals = input.getvalues();
    
    // Create the output. There are numtimeevals harmonics + the sin0 entry at the begining.
    std::vector<std::vector<densematrix>> output(numtimeevals + 1, std::vector<densematrix> {});

    // Loop on all output harmonics:
    for (long int h = 0; h < numtimeevals; h++)
    {
        // Our harmonic number is at +1 because of the sin0 term.
        long int harm = h+1;

        // The current harmonic has a frequency currentfreq*f0.
        long int currentfreq = harmonic::getfrequency(harm);
        
        // Initialise to all zero:
        densematrix currentmat(mym, myn, 0);
        double* currentvals = currentmat.getvalues();
        
        // Loop on every time step in the input matrix:
        for (long int i = 0; i < numtimeevals; i++)
        {
        	// Real part then imaginary part.
        	double coef;
        	if (harm%2 == 1)
        		coef = std::cos(2.0*pi*currentfreq*i/numtimeevals) / numtimeevals;
        	else
        		coef = std::sin(2.0*pi*currentfreq*i/numtimeevals) / numtimeevals;
        	// Correct the missing factor 2 for everything but the constant:
        	if (currentfreq > 0)
        		coef *= 2;
        
        	for (long int j = 0; j < numtransforms; j++)
        		currentvals[j] += inputvals[i*numtransforms+j] * coef;
        }

		output[harm] = {currentmat};
    }

    removeroundoffnoise(output);

	return output;
}

void myfft::removeroundoffnoise(std::vector<std::vector<densematrix>>& input, double threshold)
{    
    // First compute the max(abs()) of all harmonics:
    std::vector<double> maxabs(input.size(), 0);
    
    for (long int harm = 1; harm < input.size(); harm++)
    {
        if (input[harm].size() > 0)
            maxabs[harm] = input[harm][0].maxabs();
    }
    
    // Compute the overall harm max:
    double harmmax = 0;
    for (long int i = 1; i < input.size(); i++)
    {
        if (maxabs[i] > harmmax)
            harmmax = maxabs[i];
    }
    
    // Kill the too small harmonics:
    for (long int harm = 1; harm < input.size(); harm++)
    {
        if (input[harm].size() > 0 && (maxabs[harm] == 0 || maxabs[harm] < threshold*harmmax))
            input[harm] = {};
    }
}

densematrix myfft::inversefft(std::vector<std::vector<densematrix>>& input, long int numtimevals, long int mym, long int myn)
{
    // Get all times at which to evaluate:
    std::vector<double> timevals(numtimevals);
    
    double pi = 3.141592653589793238;
    double phasestep = 2.0*pi / ((double)(numtimevals));

    for(long int i = 0; i < numtimevals; i++)
        timevals[i] = i*phasestep;
        
    // The end result goes here. Initial value is 0.
    densematrix output(numtimevals, mym*myn, 0);
    
    // Loop on all non zero harmonics:
    densematrix sincoseval(numtimevals,1);
    double* valvec = sincoseval.getvalues();
    
    for (long int harm = 1; harm < input.size(); harm++)
    {
        if (input[harm].size() != 0)
        {
            // The current harmonic has a frequency currentfreq*f0.
            long int currentfreq = harmonic::getfrequency(harm);
        
            // Evaluate the current sin or cos term for every time step:
            if (harmonic::iscosine(harm))
            {
                for (long int i = 0; i < numtimevals; i++)
                    valvec[i] = std::cos(currentfreq*phasestep*i);
            }
            else
            {
                for (long int i = 0; i < numtimevals; i++)
                    valvec[i] = std::sin(currentfreq*phasestep*i);
            }            
            output.add( sincoseval.multiply( input[harm][0].flatten() ) );
        }
    }
    return output;
}

densematrix myfft::toelementrowformat(densematrix timestepsinrows, long int numberofelements)
{
    long int numberoftimesteps = timestepsinrows.countrows();
    long int numberofevaluationpoints = timestepsinrows.countcolumns()/numberofelements;
    densematrix output(numberofelements, numberoftimesteps*numberofevaluationpoints);
    
    double* in = timestepsinrows.getvalues();
    double* out = output.getvalues();
    
    for (long int elem = 0; elem < numberofelements; elem++)
    {    
        for (long int t = 0; t < numberoftimesteps; t++)
        {
            for (long int evalpt = 0; evalpt < numberofevaluationpoints; evalpt++)
                out[elem*numberoftimesteps*numberofevaluationpoints+t*numberofevaluationpoints+evalpt] = in[t*numberofelements*numberofevaluationpoints+elem*numberofevaluationpoints+evalpt];
        }
    }
    return output;
}

void myfft::sameharmonics(std::vector<std::vector<std::vector<densematrix>>>& notsame)
{
	if (notsame.size() <= 1)
		return;

	long int maxlen = notsame[0].size();
	for (long int i = 1; i < notsame.size(); i++)
	{
		if (notsame[i].size() > maxlen)
			maxlen = notsame[i].size();
	}
	// Make sure the length of all notsame[i] is the same:
	for (long int i = 0; i < notsame.size(); i++)
		notsame[i].resize(maxlen);
		
	for (long int h = 0; h < maxlen; h++)
	{
		long int curnumrows = -1, curnumcols = -1;
		for (long int i = 0; i < notsame.size(); i++)
		{
			if (notsame[i][h].size() == 1)
			{
				curnumrows = notsame[i][h][0].countrows();
				curnumcols = notsame[i][h][0].countcolumns();
				break;
			}
		}
		if (curnumrows == -1 || curnumcols == -1)
			continue;
			
		// Fill with full zero if an entry is missing:
		for (long int i = 0; i < notsame.size(); i++)
		{
			if (notsame[i][h].size() == 0)
				notsame[i][h] = {densematrix(curnumrows, curnumcols, 0.0)};
		}
	}
}


