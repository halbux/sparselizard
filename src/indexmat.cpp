#include "indexmat.h"
#include "gentools.h"


void indexmat::errorifempty(void)
{
    if (numrows*numcols == 0)
    {
        logs log;
        log.msg() << "Error in 'indexmat' object: cannot perform operation on empty matrix" << std::endl;
        log.error();
    }
}

indexmat::indexmat(long long int numberofrows, long long int numberofcolumns)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    myvalues = std::shared_ptr<int>(new int[numcols*numrows]);
}

indexmat::indexmat(long long int numberofrows, long long int numberofcolumns, int initvalue)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    int* myvaluesptr = new int[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = initvalue;
    
    myvalues = std::shared_ptr<int>(myvaluesptr);
}

indexmat::indexmat(long long int numberofrows, long long int numberofcolumns, std::vector<int> valvec)
{
    if (numberofrows*numberofcolumns != valvec.size())
    {
        logs log;
        log.msg() << "Error in 'indexmat' object: expected a value vector of size " << numberofrows*numberofcolumns << std::endl;
        log.error();
    }
    
    numrows = numberofrows;
    numcols = numberofcolumns;
    int* myvaluesptr = new int[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = valvec[i];
    
    myvalues = std::shared_ptr<int>(myvaluesptr);
}

indexmat::indexmat(long long int numberofrows, long long int numberofcolumns, int init, int step)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    int* myvaluesptr = new int[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = init+i*step;
    
    myvalues = std::shared_ptr<int>(myvaluesptr);
}

indexmat::indexmat(std::vector<indexmat> input)
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
            logs log;
            log.msg() << "Error in 'indexmat' object: dimension mismatch in concatenation" << std::endl;
            log.error();
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

indexmat indexmat::getresized(long long int m, long long int n)
{
    indexmat out = *this; 
    out.numrows = m;
    out.numcols = n;
    return out;
}

long long int indexmat::countpositive(void)
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

long long int indexmat::countoccurences(long long int value)
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

indexmat indexmat::removevalue(long long int toremove)
{
    int* myvaluesptr = myvalues.get();
    
    long long int num = countoccurences(toremove);
    
    indexmat output(numcols*numrows-num,1);
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

std::vector<int> indexmat::countalloccurences(int maxintval)
{
    int* myvaluesptr = myvalues.get();
    
    std::vector<int> output(maxintval+1, 0);
    
    for (long long int i = 0; i < numcols*numrows; i++)
        output[myvaluesptr[i]]++;
    
    return output;
}

std::vector<std::vector<int>> indexmat::findalloccurences(int maxintval)
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

long long int indexmat::sum(void)
{
    int* myvaluesptr = myvalues.get();
    
    long long int summed = 0;
    for (long long int i = 0; i < numcols*numrows; i++)
        summed += myvaluesptr[i];
    return summed;
}

std::vector<int> indexmat::minmax(void)
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

int indexmat::max(void)
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

void indexmat::print(void)
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

void indexmat::printsize(void)
{
    std::cout << "Matrix size is " << numrows << "x" << numcols << std::endl;
}

int* indexmat::getvalues(void)
{
    return myvalues.get();
}

indexmat indexmat::copy(void)
{
    indexmat indexmatcopy = *this;
    
    // The pointed value has to be copied as well.
    if (indexmatcopy.myvalues != NULL)
    {
        indexmatcopy.myvalues = std::shared_ptr<int>(new int[numcols*numrows]);
        int* copiedmyvaluesptr = indexmatcopy.myvalues.get();
        int* myvaluesptr = myvalues.get();
        for (long long int i = 0; i < numcols*numrows; i++)
            copiedmyvaluesptr[i] = myvaluesptr[i];
    }
        
    return indexmatcopy;
}

indexmat indexmat::gettranspose(void)
{
    indexmat output(numcols, numrows);

    int* myvaluesptr = myvalues.get();
    int* outvaluesptr = output.myvalues.get();

    for (long long int row = 0; row < numrows; row++)
    {
        for (long long int col = 0; col < numcols; col++)
            outvaluesptr[col*numrows+row] = myvaluesptr[row*numcols+col];
    }

    return output;
}

indexmat indexmat::duplicateallrowstogether(int n)
{
    indexmat output(numrows*n, numcols);
    
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

indexmat indexmat::duplicaterowsonebyone(int n)
{
    indexmat output(numrows*n, numcols);
    
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

indexmat indexmat::duplicateallcolstogether(int n)
{
    indexmat output(numrows, numcols*n);

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

indexmat indexmat::duplicatecolsonebyone(int n)
{
    indexmat output(numrows, numcols*n);

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

indexmat indexmat::extractrows(std::vector<int>& selected)
{
    long long int numselected = selected.size();

    indexmat output(numselected, numcols);

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

indexmat indexmat::extractcols(std::vector<int>& selected)
{
    long long int numselected = selected.size();

    indexmat output(numrows, numselected);

    int* vals = getvalues();
    int* outvals = output.getvalues();
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < numselected; j++)
            outvals[i*numselected+j] = vals[i*numcols+selected[j]];
    }

    return output;
}
    
indexmat indexmat::select(std::vector<bool>& sel, bool selectif)
{
    int numtrue = gentools::counttrue(sel);
    int num = numtrue;
    if (selectif == false)
        num = sel.size() - numtrue;

    indexmat output(num, 1);

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

