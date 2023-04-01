#include "densemat.h"
#include "cblas.h"


void densemat::errorifempty(void)
{
    if (numrows*numcols == 0)
    {
        logs log;
        log.msg() << "Error in 'densemat' object: cannot perform operation on empty matrix" << std::endl;
        log.error();
    }
}

densemat::densemat(long long int numberofrows, long long int numberofcolumns)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    myvalues = std::shared_ptr<double>(new double[numcols*numrows]);
}

densemat::densemat(long long int numberofrows, long long int numberofcolumns, double initvalue)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    double* myvaluesptr = new double[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = initvalue;
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

densemat::densemat(long long int numberofrows, long long int numberofcolumns, std::vector<double> valvec)
{
    if (numberofrows*numberofcolumns != valvec.size())
    {
        logs log;
        log.msg() << "Error in 'densemat' object: expected a value vector of size " << numberofrows*numberofcolumns << std::endl;
        log.error();
    }
    
    numrows = numberofrows;
    numcols = numberofcolumns;
    double* myvaluesptr = new double[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = valvec[i];
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

densemat::densemat(long long int numberofrows, long long int numberofcolumns, double init, double step)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    double* myvaluesptr = new double[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = init+i*step;
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

densemat::densemat(std::vector<densemat> input)
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
            log.msg() << "Error in 'densemat' object: dimension mismatch in concatenation" << std::endl;
            log.error();
        }
    }
    double* myvaluesptr = new double[numrows*numcols];
    
    long long int index = 0;
    for (long long int i = 0; i < input.size(); i++)
    {
        double* curvals = input[i].getvalues();
        for (long long int j = 0; j < input[i].count(); j++)
        {
            myvaluesptr[index] = curvals[j];
            index++;
        }
    }
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

void densemat::setrow(long long int rownumber, std::vector<double> rowvals)
{
    double* myvaluesptr = myvalues.get();
    for (long long int i = 0; i < numcols; i++)
        myvaluesptr[rownumber*numcols+i] = rowvals[i];
}

densemat densemat::getresized(long long int m, long long int n)
{
    densemat out = *this; 
    out.numrows = m;
    out.numcols = n;
    return out;
}

densemat densemat::getflattened(void) 
{ 
    return getresized(1, numrows*numcols);
}

void densemat::insert(long long int row, long long int col, densemat toinsert)
{
    double* myvaluesptr = myvalues.get();
    double* toinsertvaluesptr = toinsert.myvalues.get();

    for (long long int i = 0; i < toinsert.numrows; i++)
    {
        for (long long int j = 0; j < toinsert.numcols; j++)
            myvaluesptr[(row+i)*numcols+(col+j)] = toinsertvaluesptr[i*toinsert.numcols+j];
    }
}

void densemat::insertatrows(std::vector<int> selectedrows, densemat toinsert)
{
    double* myvaluesptr = myvalues.get();
    double* toinsertvaluesptr = toinsert.myvalues.get();
    
    for (long long int i = 0; i < selectedrows.size(); i++)
    {
        for (long long int col = 0; col < numcols; col++)
            myvaluesptr[selectedrows[i]*numcols+col] = toinsertvaluesptr[i*numcols+col];
    }
}

void densemat::insertatcolumns(std::vector<int> selectedcolumns, densemat toinsert)
{
    double* myvaluesptr = myvalues.get();
    double* toinsertvaluesptr = toinsert.myvalues.get();
    
    for (long long int row = 0; row < numrows; row++)
    {
        for (long long int j = 0; j < selectedcolumns.size(); j++)
            myvaluesptr[row*numcols+selectedcolumns[j]] = toinsertvaluesptr[row*toinsert.numcols+j];
    }
}

double densemat::getvalue(long long int rownumber, long long int columnnumber)
{
    double* myvaluesptr = myvalues.get();
    return myvaluesptr[rownumber*numcols+columnnumber];
}

void densemat::setvalue(long long int rownumber, long long int columnnumber, double val)
{
    double* myvaluesptr = myvalues.get();
    myvaluesptr[rownumber*numcols+columnnumber] = val;
}

void densemat::getvalues(std::vector<double>& topopulate)
{
    topopulate.resize(count());

    double* myvaluesptr = myvalues.get();
    for (long long int i = 0; i < count(); i++)
        topopulate[i] = myvaluesptr[i];
}

densemat densemat::copy(void)
{
    densemat densematcopy = *this;
    
    // The pointed value has to be copied as well.
    if (densematcopy.myvalues != NULL)
    {
        densematcopy.myvalues = std::shared_ptr<double>(new double[numcols*numrows]);
        double* copiedmyvaluesptr = densematcopy.myvalues.get();
        double* myvaluesptr = myvalues.get();
        for (long long int i = 0; i < numcols*numrows; i++)
            copiedmyvaluesptr[i] = myvaluesptr[i];
    }
        
    return densematcopy;
}

void densemat::print(void)
{
    printsize();
    
    double* myvaluesptr = myvalues.get();
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < numcols; j++)
            std::cout << myvaluesptr[i*numcols+j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void densemat::printsize(void)
{
    std::cout << "Matrix size is " << numrows << "x" << numcols << std::endl;
}

bool densemat::isallzero(void)
{
    double* myvaluesptr = myvalues.get();
    for (long long int i = 0; i < numrows*numcols; i++)
    {
        if (myvaluesptr[i] != 0)
            return false;
    }
    return true;
}

void densemat::transpose(void)
{
    istransposed = not(istransposed);
}

densemat densemat::multiply(densemat B)
{
    // Number of rows and columns of A and B after transposition:
    long long int numrowsA, numcolsA, numrowsB, numcolsB;
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
        logs log;
        log.msg() << "Error in 'densemat' object: trying to multiply a " << numrowsA << "x" << numcolsA << " matrix by a "  << numrowsB << "x" << numcolsB << std::endl;
        log.error();
    }
    
    densemat C(numrowsA, numcolsB);
    
    
    // For too small matrices the overhead of BLAS is too visible and we use a homemade product.
    long long int sizethreshold = 100*100;
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
        for (long long int m = 0; m < numrowsA; m++)
        {
            for (long long int n = 0; n < numcolsB; n++)
            {
                Cmyvaluesptr[m*numcolsB+n] = 0;
                for (long long int k = 0; k < numcolsA; k++)
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

double* densemat::getvalues(void) { return myvalues.get(); }

void densemat::addproduct(double coef, densemat B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    if (coef == 1)
    {
        for (long long int i = 0; i < numrows*numcols; i++)
            myvaluesptr[i] += Bmyvaluesptr[i];
    }
    else
    {
        for (long long int i = 0; i < numrows*numcols; i++)
            myvaluesptr[i] += coef*Bmyvaluesptr[i];
    }
}

void densemat::addproduct(densemat A, densemat B)
{
    double* myvaluesptr = myvalues.get();
    double* Amyvaluesptr = A.myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] += Amyvaluesptr[i]*Bmyvaluesptr[i];
}

densemat densemat::getproduct(double coef)
{
    densemat output(numrows, numcols);
    double* myvaluesptr = myvalues.get();
    double* outputmyvaluesptr = output.myvalues.get();
    
    if (coef == 1)
    {
        for (long long int i = 0; i < numrows*numcols; i++)
            outputmyvaluesptr[i] = myvaluesptr[i];
    }
    else
    {
        for (long long int i = 0; i < numrows*numcols; i++)
            outputmyvaluesptr[i] = coef*myvaluesptr[i];
    }
    return output;
}

densemat densemat::gettranspose(void)
{
    densemat output(numcols, numrows);
    
    double* myvaluesptr = myvalues.get();
    double* outputmyvaluesptr = output.myvalues.get();
    
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < numcols; j++)
            outputmyvaluesptr[j*numrows+i] = myvaluesptr[i*numcols+j];
    }
    
    return output;
}

densemat densemat::getinverse(void)
{
    int n = numrows;
    if (numrows != numcols)
    {
        logs log;
        log.msg() << "Error in 'densemat' object: can only invert a square densemat" << std::endl;
        log.error();
    }
    
    if (n == 0)
        return densemat(0,0);

    indexmat csrrows(1, n+1, 0,n);
    indexmat csrcols(n,n);
    int* captr = csrcols.getvalues();
    for (int r = 0; r < n; r++)
    {
        for (int c = 0; c < n; c++)
            captr[r*n+c] = c;
    }
        
    Mat bdmat;
    MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, n, n, csrrows.getvalues(), csrcols.getvalues(), myvalues.get(), &bdmat);

    MatAssemblyBegin(bdmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(bdmat, MAT_FINAL_ASSEMBLY);

    indexmat blocksizes(1,1, n);
    densemat output(n,n);
    MatInvertVariableBlockDiagonal(bdmat, 1, blocksizes.getvalues(), output.getvalues()); // output is column-major after the call

    MatDestroy(&bdmat);

    // Make row-major again:
    return output.gettranspose();
}

void densemat::multiplyelementwise(densemat B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] *= Bmyvaluesptr[i];
}

void densemat::multiplyelementwise(double val)
{
    double* myvaluesptr = myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] *= val;
}

void densemat::add(densemat B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] += Bmyvaluesptr[i];
}

void densemat::subtract(densemat B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] -= Bmyvaluesptr[i];
}

void densemat::minus(void)
{
    double* myvaluesptr = myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = -myvaluesptr[i];
}

void densemat::power(densemat exponent)
{
    double* myvaluesptr = myvalues.get();
    double* expmyvaluesptr = exponent.myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
            myvaluesptr[i] = std::pow(myvaluesptr[i], expmyvaluesptr[i]);
}

void densemat::invert(void)
{
    double* myvaluesptr = myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = 1.0/myvaluesptr[i];
}

void densemat::abs(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::abs(myvaluesptr[i]);
}

void densemat::sin(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::sin(myvaluesptr[i]);
}

void densemat::cos(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::cos(myvaluesptr[i]);
}

void densemat::tan(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::tan(myvaluesptr[i]);
}

void densemat::asin(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::asin(myvaluesptr[i]);
}

void densemat::acos(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::acos(myvaluesptr[i]);
}

void densemat::atan(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::atan(myvaluesptr[i]);
}

void densemat::log(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::log(myvaluesptr[i]);
}

void densemat::mod(double modval)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::fmod(myvaluesptr[i], modval);
}

std::vector<double> densemat::minmax(void)
{
    errorifempty();

    double* myvaluesptr = myvalues.get();
    
    double minval = myvaluesptr[0];
    double maxval = myvaluesptr[0];

    for (long long int i = 1; i < numrows*numcols; i++)
    {
        if (myvaluesptr[i] > maxval)
            maxval = myvaluesptr[i];
        if (myvaluesptr[i] < minval)
            minval = myvaluesptr[i];
    }
    return {minval, maxval};
}

double densemat::max(void)
{
    errorifempty();

    double* myvaluesptr = myvalues.get();
    
    double val = myvaluesptr[0];

    for (long long int i = 1; i < numrows*numcols; i++)
    {
        if (myvaluesptr[i] > val)
            val = myvaluesptr[i];
    }
    return val;
}

double densemat::maxabs(void)
{
    errorifempty();

    double* myvaluesptr = myvalues.get();
    
    double val = std::abs(myvaluesptr[0]);

    for (long long int i = 1; i < numrows*numcols; i++)
    {
        if (std::abs(myvaluesptr[i]) > val)
            val = std::abs(myvaluesptr[i]);
    }
    return val;
}

double densemat::sum(void)
{
    double* myvaluesptr = myvalues.get();
    double val = 0;

    for (long long int i = 0; i < numrows*numcols; i++)
            val += myvaluesptr[i];
    return val;
}

void densemat::multiplycolumns(std::vector<double> input)
{
    double* myvaluesptr = myvalues.get();
    
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < numcols; j++)
            myvaluesptr[i*numcols + j] *= input[j];
    }
}

densemat densemat::multiplyallrows(densemat input)
{
    densemat output(numrows*input.numrows, numcols);
    
    double* myvaluesptr = myvalues.get();
    double* inmyvaluesptr = input.myvalues.get();
    double* outmyvaluesptr = output.myvalues.get();
    
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < input.numrows; j++)
        {
            for (long long int k = 0; k < numcols; k++)
                outmyvaluesptr[(i*input.numrows + j)*numcols + k] = myvaluesptr[i*numcols+k] * inmyvaluesptr[j*input.numcols+k];
        }
    }
    return output;
}

densemat densemat::dofinterpoltimestf(densemat tfval)
{
    long long int fft = tfval.countrows();
    long long int gp = tfval.countcolumns();
    long long int elem = numrows;
    long long int ffd = numcols/gp;
    
    densemat output(elem, gp*ffd*fft);
    
    double* myvaluesptr = myvalues.get();
    double* tfvaluesptr = tfval.myvalues.get();
    double* outmyvaluesptr = output.myvalues.get();

    long long int ind = 0;
    for (long long int ielem = 0; ielem < elem; ielem++)
    {
        for (long long int ifft = 0; ifft < fft; ifft++)
        {
            for (long long int iffd = 0; iffd < ffd; iffd++)
            {
                for (long long int igp = 0; igp < gp; igp++)
                {
                    outmyvaluesptr[ind] = myvaluesptr[ ielem*numcols+iffd*gp+igp ] * tfvaluesptr[ifft*gp+igp];
                    ind++;
                }
            }
        }
    }

    return output;
}

void densemat::multiplycolumns(densemat input)
{
    long long int collen = input.countcolumns();
    long long int numblocks = numcols/collen;

    double* myvaluesptr = myvalues.get();
    double* invaluesptr = input.myvalues.get();
    
    long long int ind = 0;
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int b = 0; b < numblocks; b++)
        {
            for (long long int c = 0; c < collen; c++)
            {
                myvaluesptr[ind] *= invaluesptr[i*collen+c];
                ind++;
            }
        }
    }
}

densemat densemat::duplicatevertically(int n)
{
    densemat output(numrows*n, numcols);
    
    double* myvaluesptr = myvalues.get();
    double* outmyvaluesptr = output.myvalues.get();
    
    for (int duplicate = 0; duplicate < n; duplicate++)
    {
        long long int offset = duplicate*numrows*numcols;
        for (long long int i = 0; i < numrows*numcols; i++)
            outmyvaluesptr[offset + i] = myvaluesptr[i];
    }
    
    return output;
}

densemat densemat::duplicatehorizontally(int n)
{
    densemat output(numrows, numcols*n);
    
    double* myvaluesptr = myvalues.get();
    double* outmyvaluesptr = output.myvalues.get();
    
    long long int ind = 0;
    for (long long int i = 0; i < numrows; i++)
    {
        for (int duplicate = 0; duplicate < n; duplicate++)
        {
            for (long long int j = 0; j < numcols; j++)
            {
                outmyvaluesptr[ind] = myvaluesptr[i*numcols+j];
                ind++;   
            }
        }
    }
    
    return output;
}

densemat densemat::extractrows(std::vector<int>& selected)
{
    long long int numselected = selected.size();

    densemat output(numselected, numcols);
    
    double* vals = getvalues();
    double* outvals = output.getvalues();
    for (long long int i = 0; i < numselected; i++)
    {
        long long int row = selected[i];
        for (long long int j = 0; j < numcols; j++)
            outvals[i*numcols+j] = vals[row*numcols+j];
    }
    
    return output;
}

densemat densemat::extractcols(std::vector<int>& selected)
{
    long long int numselected = selected.size();
    
    densemat output(numrows, numselected);
    
    double* vals = getvalues();
    double* outvals = output.getvalues();
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < numselected; j++)
        {
            long long int col = selected[j];
            outvals[i*numselected+j] = vals[i*numcols+col];
        }
    }
    
    return output;
}

densemat densemat::extractrows(long long int rangebegin, long long int rangeend)
{
    std::vector<int> selected(rangeend-rangebegin+1);
    for (long long int i = 0; i < selected.size(); i++)
        selected[i] = rangebegin + i;
        
    return extractrows(selected);
}

densemat densemat::extractcols(long long int rangebegin, long long int rangeend)
{
    std::vector<int> selected(rangeend-rangebegin+1);
    for (long long int i = 0; i < selected.size(); i++)
        selected[i] = rangebegin + i;
        
    return extractcols(selected);
}

densemat densemat::blockdiagonaltimesvector(indexmat blocklens, densemat v)
{
    long long int nb = blocklens.count();
    
    int* blvals = blocklens.getvalues();
    double* bdvals = getvalues();
    double* vvals = v.getvalues();
    
    densemat output(v.countrows(), v.countcolumns(), 0.0);
    double* outvals = output.getvalues();
    
    // Loop on all blocks:
    long long int vecoffset = 0;
    long long int matoffset = 0;
    for (long long int n = 0; n < nb; n++)
    {
        long long int cs = blvals[n];
        // Loop on all block columns:
        for (long long int j = 0; j < cs; j++)
        {
            double curvecval = vvals[vecoffset+j];
            // Loop on all block rows:
            for (long long int i = 0; i < cs; i++)
                outvals[vecoffset+i] += bdvals[matoffset+i] * curvecval;
            matoffset += cs;
        }
        vecoffset += cs;
    }
    
    return output;
}

