#include "intdensematrix.h"


intdensematrix::intdensematrix(int numberofrows, int numberofcolumns)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    myvalues = std::shared_ptr<int>(new int[numcols*numrows]);
}

intdensematrix::intdensematrix(int numberofrows, int numberofcolumns, int initvalue)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    int* myvaluesptr = new int[numcols*numrows];
    
    for (int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = initvalue;
    
    myvalues = std::shared_ptr<int>(myvaluesptr);
}

intdensematrix::intdensematrix(int numberofrows, int numberofcolumns, int init, int step)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    int* myvaluesptr = new int[numcols*numrows];
    
    for (int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = init+i*step;
    
    myvalues = std::shared_ptr<int>(myvaluesptr);
}

int intdensematrix::countpositive(void)
{
    int* myvaluesptr = myvalues.get();
    
    int numpositive = 0;
    for (int i = 0; i < numrows*numcols; i++)
    {
        if (myvaluesptr[i] >= 0)
            numpositive++;            
    }
    return numpositive;
}

void intdensematrix::print(void)
{
    int* myvaluesptr = myvalues.get();
    
    for (int i = 0; i < numrows; i++)
    {
        for (int j = 0; j < numcols; j++)
            std::cout << myvaluesptr[i*numcols+j] << " ";
        std::cout << std::endl;
    }
}

void intdensematrix::printsize(void)
{
    std::cout << "Matrix size is " << numrows << "x" << numcols << std::endl;
}

int* intdensematrix::getvalues(void)
{
    return myvalues.get();
}

intdensematrix intdensematrix::duplicateallrowstogether(int n)
{
    intdensematrix output(numrows*n, numcols);
    
    int* myvaluesptr = myvalues.get();
    int* outvaluesptr = output.myvalues.get();
        
    int matsize = numrows*numcols;
    
    for (int duplicate = 0; duplicate < n; duplicate++)
    {
        for (int i = 0; i < matsize; i++)
            outvaluesptr[matsize*duplicate + i] = myvaluesptr[i];
    }
    
    return output;
}

intdensematrix intdensematrix::duplicaterowsonebyone(int n)
{
    intdensematrix output(numrows*n, numcols);
    
    int* myvaluesptr = myvalues.get();
    int* outvaluesptr = output.myvalues.get();
    
    for (int row = 0; row < numrows; row++)
    {
        for (int duplicate = 0; duplicate < n; duplicate++)
        {
            for (int col = 0; col < numcols; col++)
                outvaluesptr[row*(numcols*n) + duplicate*numcols + col] = myvaluesptr[row*numcols+col];
        }
    }
    
    return output;
}
    

