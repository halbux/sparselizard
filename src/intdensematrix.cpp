#include "intdensematrix.h"


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

intdensematrix::intdensematrix(long long int numberofrows, long long int numberofcolumns, const std::vector<int> valvec)
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

intdensematrix intdensematrix::transpose(void)
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

intdensematrix intdensematrix::extractrows(std::vector<int> selected)
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
    

