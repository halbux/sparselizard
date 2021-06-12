#include "intdensematrix.h"
#include "myalgorithm.h"


void intdensematrix::errorifempty(void)
{
    if (numrows*numcols == 0)
    {
        std::cout << "Error in 'intdensematrix' object: cannot perform operation on empty matrix" << std::endl;
        abort();
    }
}

intdensematrix::intdensematrix(long long int numberofrows, long long int numberofcolumns)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    myvalues = std::shared_ptr<int>(new int[numcols*numrows]);
}

intdensematrix::intdensematrix(long long int numberofrows, long long int numberofcolumns, int initvalue)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    int* myvaluesptr = new int[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = initvalue;
    
    myvalues = std::shared_ptr<int>(myvaluesptr);
}

intdensematrix::intdensematrix(long long int numberofrows, long long int numberofcolumns, std::vector<int> valvec)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    int* myvaluesptr = new int[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = valvec[i];
    
    myvalues = std::shared_ptr<int>(myvaluesptr);
}

intdensematrix::intdensematrix(long long int numberofrows, long long int numberofcolumns, int init, int step)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    int* myvaluesptr = new int[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = init+i*step;
    
    myvalues = std::shared_ptr<int>(myvaluesptr);
}

intdensematrix::intdensematrix(std::vector<intdensematrix> input)
{
    if (input.size() == 0)
        return;
        
    // Calculate the final concatenated size:
    numcols = input[0].countcolumns();
    numrows = input[0].countrows();
    for (long long int i = 1; i < input.size(); i++)
    {
        numrows += input[i].countrows();
        if (input[i].countcolumns() != numcols)
        {
            std::cout << "Error in 'intdensematrix' object: dimension mismatch in concatenation" << std::endl;
            abort();
        }
    }
    int* myvaluesptr = new int[numrows*numcols];
    
    long long int index = 0;
    for (long long int i = 0; i < input.size(); i++)
    {
        int* curvals = input[i].getvalues();
        for (long long int j = 0; j < input[i].count(); j++)
        {
            myvaluesptr[index] = curvals[j];
            index++;
        }
    }
    
    myvalues = std::shared_ptr<int>(myvaluesptr);
}

intdensematrix intdensematrix::getresized(long long int m, long long int n)
{
    intdensematrix out = *this; 
    out.numrows = m;
    out.numcols = n;
    return out;
}

long long int intdensematrix::countpositive(void)
{
    int* myvaluesptr = myvalues.get();
    
    long long int numpositive = 0;
    for (long long int i = 0; i < numrows*numcols; i++)
    {
        if (myvaluesptr[i] >= 0)
            numpositive++;            
    }
    return numpositive;
}

long long int intdensematrix::countoccurences(long long int value)
{
    int* myvaluesptr = myvalues.get();
    
    long long int num = 0;
    for (long long int i = 0; i < numcols*numrows; i++)
    {
        if (myvaluesptr[i] == value)
            num++;
    }
    return num;
}

intdensematrix intdensematrix::removevalue(long long int toremove)
{
    int* myvaluesptr = myvalues.get();
    
    long long int num = countoccurences(toremove);
    
    intdensematrix output(numcols*numrows-num,1);
    int* outvals = output.getvalues();
    
    long long int index = 0;
    for (long long int i = 0; i < numcols*numrows; i++)
    {
        if (myvaluesptr[i] != toremove)
        {
            outvals[index] = myvaluesptr[i];
            index++;
        }
    }
    return output;
}

std::vector<int> intdensematrix::countalloccurences(int maxintval)
{
    int* myvaluesptr = myvalues.get();
    
    std::vector<int> output(maxintval+1, 0);
    
    for (long long int i = 0; i < numcols*numrows; i++)
        output[myvaluesptr[i]]++;
    
    return output;
}

std::vector<std::vector<int>> intdensematrix::findalloccurences(int maxintval)
{
    int* myvaluesptr = myvalues.get();
 
    std::vector<int> numoccurences = countalloccurences(maxintval);
    
    std::vector<std::vector<int>> output(numoccurences.size());
    for (int i = 0; i < numoccurences.size(); i++)
        output[i] = std::vector<int>(numoccurences[i]);
    
    std::vector<int> indexes(numoccurences.size(), 0);
    for (long long int i = 0; i < numcols*numrows; i++)
    {
        int curvalue = myvaluesptr[i];
        output[curvalue][indexes[curvalue]] = i;
        indexes[curvalue]++;
    }
    
    return output;
}

long long int intdensematrix::sum(void)
{
    int* myvaluesptr = myvalues.get();
    
    long long int summed = 0;
    for (long long int i = 0; i < numcols*numrows; i++)
        summed += myvaluesptr[i];
    return summed;
}

std::vector<int> intdensematrix::minmax(void)
{
    errorifempty();

    int* myvaluesptr = myvalues.get();
    
    int minval = myvaluesptr[0];
    int maxval = myvaluesptr[0];

    for (long long int i = 1; i < numrows*numcols; i++)
    {
        if (myvaluesptr[i] > maxval)
            maxval = myvaluesptr[i];
        if (myvaluesptr[i] < minval)
            minval = myvaluesptr[i];
    }
    return {minval, maxval};
}

int intdensematrix::max(void)
{
    errorifempty();

    int* myvaluesptr = myvalues.get();
    
    int maxval = myvaluesptr[0];

    for (long long int i = 1; i < numrows*numcols; i++)
    {
        if (myvaluesptr[i] > maxval)
            maxval = myvaluesptr[i];
    }
    return maxval;
}

void intdensematrix::print(void)
{
    printsize();

    int* myvaluesptr = myvalues.get();
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < numcols; j++)
            std::cout << myvaluesptr[i*numcols+j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void intdensematrix::printsize(void)
{
    std::cout << "Matrix size is " << numrows << "x" << numcols << std::endl;
}

int* intdensematrix::getvalues(void)
{
    return myvalues.get();
}

intdensematrix intdensematrix::copy(void)
{
    intdensematrix intdensematrixcopy = *this;
    
    // The pointed value has to be copied as well.
    if (intdensematrixcopy.myvalues != NULL)
    {
        intdensematrixcopy.myvalues = std::shared_ptr<int>(new int[numcols*numrows]);
        int* copiedmyvaluesptr = intdensematrixcopy.myvalues.get();
        int* myvaluesptr = myvalues.get();
        for (long long int i = 0; i < numcols*numrows; i++)
            copiedmyvaluesptr[i] = myvaluesptr[i];
    }
        
    return intdensematrixcopy;
}

intdensematrix intdensematrix::gettranspose(void)
{
    intdensematrix output(numcols, numrows);

    int* myvaluesptr = myvalues.get();
    int* outvaluesptr = output.myvalues.get();

    for (long long int row = 0; row < numrows; row++)
    {
        for (long long int col = 0; col < numcols; col++)
            outvaluesptr[col*numrows+row] = myvaluesptr[row*numcols+col];
    }

    return output;
}

intdensematrix intdensematrix::duplicateallrowstogether(int n)
{
    intdensematrix output(numrows*n, numcols);
    
    int* myvaluesptr = myvalues.get();
    int* outvaluesptr = output.myvalues.get();
        
    long long int matsize = numrows*numcols;
    
    for (int duplicate = 0; duplicate < n; duplicate++)
    {
        for (long long int i = 0; i < matsize; i++)
            outvaluesptr[matsize*duplicate + i] = myvaluesptr[i];
    }
    
    return output;
}

intdensematrix intdensematrix::duplicaterowsonebyone(int n)
{
    intdensematrix output(numrows*n, numcols);
    
    int* myvaluesptr = myvalues.get();
    int* outvaluesptr = output.myvalues.get();
    
    for (long long int row = 0; row < numrows; row++)
    {
        for (int duplicate = 0; duplicate < n; duplicate++)
        {
            for (long long int col = 0; col < numcols; col++)
                outvaluesptr[row*(numcols*n) + duplicate*numcols + col] = myvaluesptr[row*numcols+col];
        }
    }
    
    return output;
}

intdensematrix intdensematrix::duplicateallcolstogether(int n)
{
    intdensematrix output(numrows, numcols*n);

    int* myvaluesptr = myvalues.get();
    int* outvaluesptr = output.myvalues.get();

    long long int ind = 0;
    for (long long int i = 0; i < numrows; i++)
    {
        for (int duplicate = 0; duplicate < n; duplicate++)
        {
            for (long long int j = 0; j < numcols; j++)
            {
                outvaluesptr[ind] = myvaluesptr[i*numcols+j];
                ind++;   
            }
        }
    }

    return output;
}

intdensematrix intdensematrix::duplicatecolsonebyone(int n)
{
    intdensematrix output(numrows, numcols*n);

    int* myvaluesptr = myvalues.get();
    int* outvaluesptr = output.myvalues.get();

    long long int ind = 0;
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < numcols; j++)
        {
            for (int duplicate = 0; duplicate < n; duplicate++)
            {
                outvaluesptr[ind] = myvaluesptr[i*numcols+j];
                ind++;   
            }
        }
    }

    return output;
}

intdensematrix intdensematrix::extractrows(std::vector<int>& selected)
{
    long long int numselected = selected.size();

    intdensematrix output(numselected, numcols);

    int* vals = getvalues();
    int* outvals = output.getvalues();
    for (long long int i = 0; i < numselected; i++)
    {
        long long int row = selected[i];
        for (long long int j = 0; j < numcols; j++)
            outvals[i*numcols+j] = vals[row*numcols+j];
    }

    return output;
}

intdensematrix intdensematrix::extractcols(std::vector<int>& selected)
{
    long long int numselected = selected.size();

    intdensematrix output(numrows, numselected);

    int* vals = getvalues();
    int* outvals = output.getvalues();
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < numselected; j++)
            outvals[i*numselected+j] = vals[i*numcols+selected[j]];
    }

    return output;
}
    
intdensematrix intdensematrix::select(std::vector<bool>& sel, bool selectif)
{
    int numtrue = myalgorithm::counttrue(sel);
    int num = numtrue;
    if (selectif == false)
        num = sel.size() - numtrue;

    intdensematrix output(num, 1);

    int* myvaluesptr = myvalues.get();
    int* outvaluesptr = output.myvalues.get();

    int index = 0;
    for (int i = 0; i < sel.size(); i++)
    {
        if (sel[i] == selectif)
        {
            outvaluesptr[index] = myvaluesptr[i];
            index++;
        }
    }

    return output;
}

