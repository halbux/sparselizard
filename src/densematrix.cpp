#include "densematrix.h"


void densematrix::errorifempty(void)
{
    if (numrows*numcols == 0)
    {
        std::cout << "Error in 'densematrix' object: cannot perform operation on empty matrix" << std::endl;
        abort();
    }
}

densematrix::densematrix(int numberofrows, int numberofcolumns)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    myvalues = std::shared_ptr<double>(new double[numcols*numrows]);
}

densematrix::densematrix(int numberofrows, int numberofcolumns, double initvalue)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    double* myvaluesptr = new double[numcols*numrows];
    
    for (int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = initvalue;
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

densematrix::densematrix(int numberofrows, int numberofcolumns, const std::vector<double> valvec)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    double* myvaluesptr = new double[numcols*numrows];
    
    for (int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = valvec[i];
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

densematrix::densematrix(int numberofrows, int numberofcolumns, double init, double step)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    double* myvaluesptr = new double[numcols*numrows];
    
    for (int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = init+i*step;
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

densematrix::densematrix(std::vector<densematrix> input)
{
    if (input.size() == 0)
        return;
        
    // Calculate the final concatenated size:
    numcols = input[0].countcolumns();
    numrows = input[0].countrows();
    for (int i = 1; i < input.size(); i++)
    {
        numrows += input[i].countrows();
        if (input[i].countcolumns() != numcols)
        {
            std::cout << "Error in 'densematrix' object: dimension mismatch in concatenation" << std::endl;
            abort();
        }
    }
    double* myvaluesptr = new double[numrows*numcols];
    
    int index = 0;
    for (int i = 0; i < input.size(); i++)
    {
        double* curvals = input[i].getvalues();
        for (int j = 0; j < input[i].count(); j++)
        {
            myvaluesptr[index] = curvals[j];
            index++;
        }
    }
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

void densematrix::setrow(int rownumber, std::vector<double> rowvals)
{
    double* myvaluesptr = myvalues.get();
    for (int i = 0; i < numcols; i++)
        myvaluesptr[rownumber*numcols+i] = rowvals[i];
}

densematrix densematrix::flatten(void) 
{ 
    densematrix out = *this; 
    out.numcols = numrows*numcols; 
    out.numrows = 1; 
    return out; 
}


void densematrix::insert(int row, int col, densematrix toinsert)
{
    double* myvaluesptr = myvalues.get();
    double* toinsertvaluesptr = toinsert.myvalues.get();

    for (int i = 0; i < toinsert.numrows; i++)
    {
        for (int j = 0; j < toinsert.numcols; j++)
            myvaluesptr[(row+i)*numcols+(col+j)] = toinsertvaluesptr[i*toinsert.numcols+j];
    }
}

void densematrix::insertatrows(std::vector<int> selectedrows, densematrix toinsert)
{
    double* myvaluesptr = myvalues.get();
    double* toinsertvaluesptr = toinsert.myvalues.get();
    
    for (int i = 0; i < selectedrows.size(); i++)
    {
        for (int col = 0; col < numcols; col++)
            myvaluesptr[selectedrows[i]*numcols+col] = toinsertvaluesptr[i*numcols+col];
    }
}

void densematrix::insertatcolumns(std::vector<int> selectedcolumns, densematrix toinsert)
{
    double* myvaluesptr = myvalues.get();
    double* toinsertvaluesptr = toinsert.myvalues.get();
    
    for (int row = 0; row < numrows; row++)
    {
        for (int j = 0; j < selectedcolumns.size(); j++)
            myvaluesptr[row*numcols+selectedcolumns[j]] = toinsertvaluesptr[row*toinsert.numcols+j];
    }
}

double densematrix::getvalue(int rownumber, int columnnumber)
{
    double* myvaluesptr = myvalues.get();
    return myvaluesptr[rownumber*numcols+columnnumber];
}

void densematrix::setvalue(int rownumber, int columnnumber, double val)
{
    double* myvaluesptr = myvalues.get();
    myvaluesptr[rownumber*numcols+columnnumber] = val;
}

void densematrix::getvalues(std::vector<double>& topopulate)
{
    topopulate.resize(count());

    double* myvaluesptr = myvalues.get();
    for (int i = 0; i < count(); i++)
        topopulate[i] = myvaluesptr[i];
}

densematrix densematrix::copy(void)
{
    densematrix densematrixcopy = *this;
    
    // The pointed value has to be copied as well.
    if (densematrixcopy.myvalues != NULL)
    {
        densematrixcopy.myvalues = std::shared_ptr<double>(new double[numcols*numrows]);
        double* copiedmyvaluesptr = densematrixcopy.myvalues.get();
        double* myvaluesptr = myvalues.get();
        for (int i = 0; i < numcols*numrows; i++)
            copiedmyvaluesptr[i] = myvaluesptr[i];
    }
        
    return densematrixcopy;
}

void densematrix::print(void)
{
    printsize();
    
    double* myvaluesptr = myvalues.get();
    for (int i = 0; i < numrows; i++)
    {
        for (int j = 0; j < numcols; j++)
            std::cout << myvaluesptr[i*numcols+j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void densematrix::printsize(void)
{
    std::cout << "Matrix size is " << numrows << "x" << numcols << std::endl;
}

bool densematrix::isallzero(void)
{
    double* myvaluesptr = myvalues.get();
    for (int i = 0; i < numrows*numcols; i++)
    {
        if (myvaluesptr[i] != 0)
            return false;
    }
    return true;
}

void densematrix::transpose(void)
{
    istransposed = not(istransposed);
}

densematrix densematrix::multiply(densematrix B)
{
    // Number of rows and columns of A and B after transposition:
    int numrowsA, numcolsA, numrowsB, numcolsB;
    if (not(istransposed))
    {
        numrowsA = numrows;
        numcolsA = numcols;
    }
    else
    {
        numrowsA = numcols;
        numcolsA = numrows;
    }
    if (not(B.istransposed))
    {
        numrowsB = B.numrows;
        numcolsB = B.numcols;
    }
    else
    {
        numrowsB = B.numcols;
        numcolsB = B.numrows;
    }
    
    if (numcolsA != numrowsB)
    {
        std::cout << "Error on 'densematrix' object: trying to multiply a " << numrowsA << "x" << numcolsA << " matrix by a "  << numrowsB << "x" << numcolsB << std::endl;
        abort();
    }
    
    densematrix C(numrowsA, numcolsB);
    
    
    // For too small matrices the overhead of BLAS is too visible and we use a homemade product.
    int sizethreshold = 100*100;
    if (numrowsA*numcolsA > sizethreshold || numrowsB*numcolsB > sizethreshold)
    {
        // 'cblas_dgemm' computes alpha*A*B+beta*C and puts the result in C. Here we only compute A*B thus:
        double alpha = 1, beta = 0;

        if (    istransposed  &&     B.istransposed)
            cblas_dgemm(CblasRowMajor, CblasTrans,   CblasTrans,   numrowsA, numcolsB, numcolsA, alpha, myvalues.get(), numcols, B.myvalues.get(), B.numcols, beta, C.myvalues.get(), numcolsB);
        if (    istransposed  && not(B.istransposed))
            cblas_dgemm(CblasRowMajor, CblasTrans,   CblasNoTrans, numrowsA, numcolsB, numcolsA, alpha, myvalues.get(), numcols, B.myvalues.get(), B.numcols, beta, C.myvalues.get(), numcolsB);
        if (not(istransposed) &&     B.istransposed)
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,   numrowsA, numcolsB, numcolsA, alpha, myvalues.get(), numcols, B.myvalues.get(), B.numcols, beta, C.myvalues.get(), numcolsB);
        if (not(istransposed) && not(B.istransposed))
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numrowsA, numcolsB, numcolsA, alpha, myvalues.get(), numcols, B.myvalues.get(), B.numcols, beta, C.myvalues.get(), numcolsB);
    }
    else
    {
        double* myvaluesptr = myvalues.get();
        double* Bmyvaluesptr = B.myvalues.get();
        double* Cmyvaluesptr = C.myvalues.get();
        
        // Hand coded matrix matrix product:
        for (int m = 0; m < numrowsA; m++)
        {
            for (int n = 0; n < numcolsB; n++)
            {
                Cmyvaluesptr[m*numcolsB+n] = 0;
                for (int k = 0; k < numcolsA; k++)
                {
                    double valA, valB;
                    if (not(istransposed))
                        valA = myvaluesptr[m*numcols+k];
                    else
                        valA = myvaluesptr[k*numcols+m];
                    
                    if (not(B.istransposed))
                        valB = Bmyvaluesptr[k*B.numcols+n];
                    else
                        valB = Bmyvaluesptr[n*B.numcols+k];
                    
                    Cmyvaluesptr[m*numcolsB+n] += valA*valB;
                }
            }
        }
    }
    
    return C;
}

double* densematrix::getvalues(void) { return myvalues.get(); }

void densematrix::addproduct(double coef, densematrix B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    if (coef == 1)
    {
        for (int i = 0; i < numrows*numcols; i++)
            myvaluesptr[i] += Bmyvaluesptr[i];
    }
    else
    {
        for (int i = 0; i < numrows*numcols; i++)
            myvaluesptr[i] += coef*Bmyvaluesptr[i];
    }
}

densematrix densematrix::returnproduct(double coef)
{
    densematrix output(numrows, numcols);
    double* myvaluesptr = myvalues.get();
    double* outputmyvaluesptr = output.myvalues.get();
    
    if (coef == 1)
    {
        for (int i = 0; i < numrows*numcols; i++)
            outputmyvaluesptr[i] = myvaluesptr[i];
    }
    else
    {
        for (int i = 0; i < numrows*numcols; i++)
            outputmyvaluesptr[i] = coef*myvaluesptr[i];
    }
    return output;
}

void densematrix::multiplyelementwise(densematrix B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] *= Bmyvaluesptr[i];
}

void densematrix::multiplyelementwise(double val)
{
    double* myvaluesptr = myvalues.get();
    
    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] *= val;
}

void densematrix::add(densematrix B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] += Bmyvaluesptr[i];
}

void densematrix::subtract(densematrix B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] -= Bmyvaluesptr[i];
}

void densematrix::minus(void)
{
    double* myvaluesptr = myvalues.get();
    
    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = -myvaluesptr[i];
}

void densematrix::power(densematrix exponent)
{
    double* myvaluesptr = myvalues.get();
    double* expmyvaluesptr = exponent.myvalues.get();
    
    for (int i = 0; i < numrows*numcols; i++)
            myvaluesptr[i] = std::pow(myvaluesptr[i], expmyvaluesptr[i]);
}

void densematrix::invert(void)
{
    double* myvaluesptr = myvalues.get();
    
    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = 1.0/myvaluesptr[i];
}

void densematrix::abs(void)
{
    double* myvaluesptr = myvalues.get();

    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::abs(myvaluesptr[i]);
}

void densematrix::sin(void)
{
    double* myvaluesptr = myvalues.get();

    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::sin(myvaluesptr[i]);
}

void densematrix::cos(void)
{
    double* myvaluesptr = myvalues.get();

    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::cos(myvaluesptr[i]);
}

void densematrix::tan(void)
{
    double* myvaluesptr = myvalues.get();

    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::tan(myvaluesptr[i]);
}

void densematrix::asin(void)
{
    double* myvaluesptr = myvalues.get();

    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::asin(myvaluesptr[i]);
}

void densematrix::acos(void)
{
    double* myvaluesptr = myvalues.get();

    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::acos(myvaluesptr[i]);
}

void densematrix::atan(void)
{
    double* myvaluesptr = myvalues.get();

    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::atan(myvaluesptr[i]);
}

void densematrix::log10(void)
{
    double* myvaluesptr = myvalues.get();

    for (int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::log10(myvaluesptr[i]);
}

double densematrix::max(void)
{
    errorifempty();

    double* myvaluesptr = myvalues.get();
    
    double val = myvaluesptr[0];

    for (int i = 1; i < numrows*numcols; i++)
    {
        if (myvaluesptr[i] > val)
            val = myvaluesptr[i];
    }
    return val;
}

double densematrix::maxabs(void)
{
    errorifempty();

    double* myvaluesptr = myvalues.get();
    
    double val = std::abs(myvaluesptr[0]);

    for (int i = 1; i < numrows*numcols; i++)
    {
        if (std::abs(myvaluesptr[i]) > val)
            val = std::abs(myvaluesptr[i]);
    }
    return val;
}

double densematrix::sum(void)
{
    errorifempty();
    
    double* myvaluesptr = myvalues.get();
    double val = 0;

    for (int i = 0; i < numrows*numcols; i++)
            val += myvaluesptr[i];
    return val;
}

void densematrix::multiplycolumns(std::vector<double> input)
{
    double* myvaluesptr = myvalues.get();
    
    for (int i = 0; i < numrows; i++)
    {
        for (int j = 0; j < numcols; j++)
            myvaluesptr[i*numcols + j] *= input[j];
    }
}

densematrix densematrix::multiplyallrows(densematrix input)
{
    densematrix output(numrows*input.numrows, numcols);
    
    double* myvaluesptr = myvalues.get();
    double* inmyvaluesptr = input.myvalues.get();
    double* outmyvaluesptr = output.myvalues.get();
    
    for (int i = 0; i < numrows; i++)
    {
        for (int j = 0; j < input.numrows; j++)
        {
            for (int k = 0; k < numcols; k++)
                outmyvaluesptr[(i*input.numrows + j)*numcols + k] = myvaluesptr[i*numcols+k] * inmyvaluesptr[j*input.numcols+k];
        }
    }
    return output;
}

densematrix densematrix::duplicatevertically(int n)
{
    densematrix output(numrows*n, numcols);
    
    double* myvaluesptr = myvalues.get();
    double* outmyvaluesptr = output.myvalues.get();
    
    for (int duplicate = 0; duplicate < n; duplicate++)
    {
        int offset = duplicate*numrows*numcols;
        for (int i = 0; i < numrows*numcols; i++)
            outmyvaluesptr[offset + i] = myvaluesptr[i];
    }
    
    return output;
}

densematrix densematrix::extractrows(std::vector<int> selected)
{
    int numselected = selected.size();

    densematrix output(numselected, numcols);
    
    double* vals = getvalues();
    double* outvals = output.getvalues();
    for (int i = 0; i < numselected; i++)
    {
        int row = selected[i];
        for (int j = 0; j < numcols; j++)
            outvals[i*numcols+j] = vals[row*numcols+j];
    }
    
    return output;
}

densematrix densematrix::extractcols(std::vector<int> selected)
{
    int numselected = selected.size();
    
    densematrix output(numrows, numselected);
    
    double* vals = getvalues();
    double* outvals = output.getvalues();
    for (int i = 0; i < numrows; i++)
    {
        for (int j = 0; j < numselected; j++)
        {
            int col = selected[j];
            outvals[i*numselected+j] = vals[i*numcols+col];
        }
    }
    
    return output;
}

densematrix densematrix::extractrows(int rangebegin, int rangeend)
{
    std::vector<int> selected(rangeend-rangebegin+1);
    for (int i = 0; i < selected.size(); i++)
        selected[i] = rangebegin + i;
        
    return extractrows(selected);
}

densematrix densematrix::extractcols(int rangebegin, int rangeend)
{
    std::vector<int> selected(rangeend-rangebegin+1);
    for (int i = 0; i < selected.size(); i++)
        selected[i] = rangebegin + i;
        
    return extractcols(selected);
}

