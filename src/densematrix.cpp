#include "densematrix.h"
#include "cblas.h"


void densematrix::errorifempty(void)
{
    if (numrows*numcols == 0)
    {
        std::cout << "Error in 'densematrix' object: cannot perform operation on empty matrix" << std::endl;
        abort();
    }
}

densematrix::densematrix(long long int numberofrows, long long int numberofcolumns)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    myvalues = std::shared_ptr<double>(new double[numcols*numrows]);
}

densematrix::densematrix(long long int numberofrows, long long int numberofcolumns, double initvalue)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    double* myvaluesptr = new double[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = initvalue;
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

densematrix::densematrix(long long int numberofrows, long long int numberofcolumns, const std::vector<double> valvec)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    double* myvaluesptr = new double[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
        myvaluesptr[i] = valvec[i];
    
    myvalues = std::shared_ptr<double>(myvaluesptr);
}

densematrix::densematrix(long long int numberofrows, long long int numberofcolumns, double init, double step)
{
    numrows = numberofrows;
    numcols = numberofcolumns;
    double* myvaluesptr = new double[numcols*numrows];
    
    for (long long int i = 0; i < numcols*numrows; i++)
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
    for (long long int i = 1; i < input.size(); i++)
    {
        numrows += input[i].countrows();
        if (input[i].countcolumns() != numcols)
        {
            std::cout << "Error in 'densematrix' object: dimension mismatch in concatenation" << std::endl;
            abort();
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

void densematrix::setrow(long long int rownumber, std::vector<double> rowvals)
{
    double* myvaluesptr = myvalues.get();
    for (long long int i = 0; i < numcols; i++)
        myvaluesptr[rownumber*numcols+i] = rowvals[i];
}

densematrix densematrix::flatten(void) 
{ 
    densematrix out = *this; 
    out.numcols = numrows*numcols; 
    out.numrows = 1; 
    return out; 
}


void densematrix::insert(long long int row, long long int col, densematrix toinsert)
{
    double* myvaluesptr = myvalues.get();
    double* toinsertvaluesptr = toinsert.myvalues.get();

    for (long long int i = 0; i < toinsert.numrows; i++)
    {
        for (long long int j = 0; j < toinsert.numcols; j++)
            myvaluesptr[(row+i)*numcols+(col+j)] = toinsertvaluesptr[i*toinsert.numcols+j];
    }
}

void densematrix::insertatrows(std::vector<int> selectedrows, densematrix toinsert)
{
    double* myvaluesptr = myvalues.get();
    double* toinsertvaluesptr = toinsert.myvalues.get();
    
    for (long long int i = 0; i < selectedrows.size(); i++)
    {
        for (long long int col = 0; col < numcols; col++)
            myvaluesptr[selectedrows[i]*numcols+col] = toinsertvaluesptr[i*numcols+col];
    }
}

void densematrix::insertatcolumns(std::vector<int> selectedcolumns, densematrix toinsert)
{
    double* myvaluesptr = myvalues.get();
    double* toinsertvaluesptr = toinsert.myvalues.get();
    
    for (long long int row = 0; row < numrows; row++)
    {
        for (long long int j = 0; j < selectedcolumns.size(); j++)
            myvaluesptr[row*numcols+selectedcolumns[j]] = toinsertvaluesptr[row*toinsert.numcols+j];
    }
}

double densematrix::getvalue(long long int rownumber, long long int columnnumber)
{
    double* myvaluesptr = myvalues.get();
    return myvaluesptr[rownumber*numcols+columnnumber];
}

void densematrix::setvalue(long long int rownumber, long long int columnnumber, double val)
{
    double* myvaluesptr = myvalues.get();
    myvaluesptr[rownumber*numcols+columnnumber] = val;
}

void densematrix::getvalues(std::vector<double>& topopulate)
{
    topopulate.resize(count());

    double* myvaluesptr = myvalues.get();
    for (long long int i = 0; i < count(); i++)
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
        for (long long int i = 0; i < numcols*numrows; i++)
            copiedmyvaluesptr[i] = myvaluesptr[i];
    }
        
    return densematrixcopy;
}

void densematrix::print(void)
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

void densematrix::printsize(void)
{
    std::cout << "Matrix size is " << numrows << "x" << numcols << std::endl;
}

bool densematrix::isallzero(void)
{
    double* myvaluesptr = myvalues.get();
    for (long long int i = 0; i < numrows*numcols; i++)
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
        std::cout << "Error on 'densematrix' object: trying to multiply a " << numrowsA << "x" << numcolsA << " matrix by a "  << numrowsB << "x" << numcolsB << std::endl;
        abort();
    }
    
    densematrix C(numrowsA, numcolsB);
    
    
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

double* densematrix::getvalues(void) { return myvalues.get(); }

void densematrix::addproduct(double coef, densematrix B)
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

void densematrix::addproduct(densematrix A, densematrix B)
{
    double* myvaluesptr = myvalues.get();
    double* Amyvaluesptr = A.myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] += Amyvaluesptr[i]*Bmyvaluesptr[i];
}

densematrix densematrix::returnproduct(double coef)
{
    densematrix output(numrows, numcols);
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

void densematrix::multiplyelementwise(densematrix B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] *= Bmyvaluesptr[i];
}

void densematrix::multiplyelementwise(double val)
{
    double* myvaluesptr = myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] *= val;
}

void densematrix::add(densematrix B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] += Bmyvaluesptr[i];
}

void densematrix::subtract(densematrix B)
{
    double* myvaluesptr = myvalues.get();
    double* Bmyvaluesptr = B.myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] -= Bmyvaluesptr[i];
}

void densematrix::minus(void)
{
    double* myvaluesptr = myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = -myvaluesptr[i];
}

void densematrix::power(densematrix exponent)
{
    double* myvaluesptr = myvalues.get();
    double* expmyvaluesptr = exponent.myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
            myvaluesptr[i] = std::pow(myvaluesptr[i], expmyvaluesptr[i]);
}

void densematrix::invert(void)
{
    double* myvaluesptr = myvalues.get();
    
    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = 1.0/myvaluesptr[i];
}

void densematrix::abs(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::abs(myvaluesptr[i]);
}

void densematrix::sin(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::sin(myvaluesptr[i]);
}

void densematrix::cos(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::cos(myvaluesptr[i]);
}

void densematrix::tan(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::tan(myvaluesptr[i]);
}

void densematrix::asin(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::asin(myvaluesptr[i]);
}

void densematrix::acos(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::acos(myvaluesptr[i]);
}

void densematrix::atan(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::atan(myvaluesptr[i]);
}

void densematrix::log10(void)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::log10(myvaluesptr[i]);
}

void densematrix::mod(double modval)
{
    double* myvaluesptr = myvalues.get();

    for (long long int i = 0; i < numrows*numcols; i++)
        myvaluesptr[i] = std::fmod(myvaluesptr[i], modval);
}

std::vector<double> densematrix::minmax(void)
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

double densematrix::max(void)
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

double densematrix::maxabs(void)
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

double densematrix::sum(void)
{
    errorifempty();
    
    double* myvaluesptr = myvalues.get();
    double val = 0;

    for (long long int i = 0; i < numrows*numcols; i++)
            val += myvaluesptr[i];
    return val;
}

void densematrix::multiplycolumns(std::vector<double> input)
{
    double* myvaluesptr = myvalues.get();
    
    for (long long int i = 0; i < numrows; i++)
    {
        for (long long int j = 0; j < numcols; j++)
            myvaluesptr[i*numcols + j] *= input[j];
    }
}

densematrix densematrix::multiplyallrows(densematrix input)
{
    densematrix output(numrows*input.numrows, numcols);
    
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

densematrix densematrix::dofinterpoltimestf(densematrix tfval)
{
    long long int fft = tfval.countrows();
    long long int gp = tfval.countcolumns();
    long long int elem = numrows;
    long long int ffd = numcols/gp;
    
    densematrix output(elem, gp*ffd*fft);
    
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

void densematrix::multiplycolumns(densematrix input)
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

densematrix densematrix::duplicatevertically(int n)
{
    densematrix output(numrows*n, numcols);
    
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

densematrix densematrix::duplicatehorizontally(int n)
{
    densematrix output(numrows, numcols*n);
    
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

densematrix densematrix::extractrows(std::vector<int> selected)
{
    long long int numselected = selected.size();

    densematrix output(numselected, numcols);
    
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

densematrix densematrix::extractcols(std::vector<int> selected)
{
    long long int numselected = selected.size();
    
    densematrix output(numrows, numselected);
    
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

densematrix densematrix::extractrows(long long int rangebegin, long long int rangeend)
{
    std::vector<int> selected(rangeend-rangebegin+1);
    for (long long int i = 0; i < selected.size(); i++)
        selected[i] = rangebegin + i;
        
    return extractrows(selected);
}

densematrix densematrix::extractcols(long long int rangebegin, long long int rangeend)
{
    std::vector<int> selected(rangeend-rangebegin+1);
    for (long long int i = 0; i < selected.size(); i++)
        selected[i] = rangebegin + i;
        
    return extractcols(selected);
}

densematrix densematrix::blockdiagonaltimesvector(intdensematrix blocklens, densematrix v)
{
    long long int nb = blocklens.count();
    
    int* blvals = blocklens.getvalues();
    double* bdvals = getvalues();
    double* vvals = v.getvalues();
    
    densematrix output(v.countrows(), v.countcolumns(), 0.0);
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

