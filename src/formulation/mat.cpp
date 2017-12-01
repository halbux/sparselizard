#include "mat.h"


void mat::errorifpointerisnull(void)
{
    if (rawmatptr == NULL)
    {
        std::cout << "Error in 'mat' object: cannot perform the operation (matrix is undefined)" << std::endl;
        abort();
    }
}

mat::mat(formulation myformulation, intdensematrix rowadresses, intdensematrix coladresses, densematrix vals)
{
	rawmatptr = shared_ptr<rawmat>(new rawmat(myformulation.getdofmanager()));
	rawmatptr->accumulate(rowadresses, coladresses, vals);
	rawmatptr->process();
}

int mat::countrows(void) { errorifpointerisnull(); return rawmatptr->countrows(); }
int mat::countcolumns(void) { errorifpointerisnull(); return rawmatptr->countcolumns(); }
        
int mat::countnnz(void) { errorifpointerisnull(); return rawmatptr->countnnz(); }

void mat::removeconstraints(void) { errorifpointerisnull(); rawmatptr->removeconstraints(); };

void mat::reuselu(void) { errorifpointerisnull(); rawmatptr->reuselu(); }

Mat mat::getpetsc(void) { errorifpointerisnull(); return rawmatptr->getpetsc(); }   

void mat::print(void) { errorifpointerisnull(); rawmatptr->print(); }

mat mat::copy(void)
{
    Mat A = getpetsc();
    Mat output;
    MatConvert(A, MATSAME, MAT_INITIAL_MATRIX, &output);
    return mat(shared_ptr<rawmat>(new rawmat(  rawmatptr->getdofmanager(), output  )));
}


mat mat::operator+(void) { return copy(); }
mat mat::operator-(void) { return *this*-1; }

mat mat::operator*(double input)
{
    Mat A = getpetsc();
    Mat output;
    MatDuplicate(A, MAT_SHARE_NONZERO_PATTERN, &output);
    MatAXPY(output, input, A, SAME_NONZERO_PATTERN);
    return mat(shared_ptr<rawmat>(new rawmat(  rawmatptr->getdofmanager(), output  )));
}

mat mat::operator*(mat input)
{
    Mat A = getpetsc();
    Mat B = input.getpetsc();
    Mat output;
    MatMatMult(A, B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &output);
    return mat(shared_ptr<rawmat>(new rawmat(  rawmatptr->getdofmanager(), output  )));
}

mat mat::operator+(mat input)
{
    mat copied = copy();
    Mat Y = copied.getpetsc();
    Mat X = input.getpetsc();
    MatAXPY(Y, 1, X, DIFFERENT_NONZERO_PATTERN);
    return copied;
}

mat mat::operator-(mat input)
{
    mat copied = copy();
    Mat Y = copied.getpetsc();
    Mat X = input.getpetsc();
    MatAXPY(Y, -1, X, DIFFERENT_NONZERO_PATTERN);
    return copied;
}

vec mat::operator*(vec input)
{
    Mat A = getpetsc();
    Vec x = input.getpetsc();
    Vec output;
    VecDuplicate(x, &output);
    MatMult(A, x, output);
    return vec(shared_ptr<rawvec>(new rawvec(  rawmatptr->getdofmanager(), output  )));
}


mat operator*(double inputdouble, mat inputmat) { return inputmat*inputdouble; }



