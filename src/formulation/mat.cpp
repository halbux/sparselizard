#include "mat.h"


void mat::errorifpointerisnull(void)
{
    if (rawmatptr == NULL)
    {
        logs log;
        log.msg() << "Error in 'mat' object: cannot perform the operation (matrix is undefined)" << std::endl;
        log.error();
    }
}

void mat::errorifinvalidated(void)
{
    if (rawmatptr != NULL && rawmatptr->getdofmanager()->ismanaged() && rawmatptr->getmeshnumber() != universe::getrawmesh()->getmeshnumber())
    {
        logs log;
        log.msg() << "Error in 'mat' object: matrix cannot be used anymore (invalidated by hp-adaptivity)" << std::endl;
        log.error();
    }
}

mat::mat(long long int matsize, indexmat rowadresses, indexmat coladresses, densemat vals)
{
    rawmatptr = std::shared_ptr<rawmat>(new rawmat(std::shared_ptr<dofmanager>(new dofmanager(matsize))));
    rawmatptr->accumulate(rowadresses, coladresses, vals);
    std::vector<bool> isconstr(matsize, false);
    rawmatptr->process(isconstr);
    rawmatptr->clearfragments();
}

mat::mat(formulation myformulation, indexmat rowadresses, indexmat coladresses, densemat vals)
{
    rawmatptr = std::shared_ptr<rawmat>(new rawmat(myformulation.getdofmanager()));
    rawmatptr->accumulate(rowadresses, coladresses, vals);
    std::vector<bool> isconstr = myformulation.getdofmanager()->isconstrained();
    rawmatptr->process(isconstr);
    rawmatptr->clearfragments();
}

long long int mat::countrows(void) { errorifpointerisnull(); errorifinvalidated(); return rawmatptr->countrows(); }
long long int mat::countcolumns(void) { errorifpointerisnull(); errorifinvalidated(); return rawmatptr->countcolumns(); }
        
long long int mat::countnnz(void) { errorifpointerisnull(); errorifinvalidated(); return rawmatptr->countnnz(); }

void mat::reusefactorization(void) { errorifpointerisnull(); errorifinvalidated(); rawmatptr->reusefactorization(); }

std::shared_ptr<rawmat> mat::getpointer(void)
{
    errorifinvalidated();
    return rawmatptr;
}
        
vec mat::xbmerge(vec x, vec b)
{
    errorifpointerisnull(); errorifinvalidated();

    vec output(std::shared_ptr<rawvec>(new rawvec(b.getpointer()->getdofmanager())));
    indexmat ainds = getainds();
    indexmat dinds = getdinds();
    densemat xvals = x.getallvalues();
    densemat bdvals = b.getvalues(dinds);
    output.setvalues(ainds, xvals);
    output.setvalues(dinds, bdvals);
    return output;
}

vec mat::x0merge(vec x)
{
    errorifpointerisnull(); errorifinvalidated();
    
    vec output(std::shared_ptr<rawvec>(new rawvec(rawmatptr->getdofmanager())));
    output.setvalues(getainds(), x.getallvalues());
    return output;
}

vec mat::eliminate(vec b)
{
    errorifpointerisnull(); errorifinvalidated();
    
    indexmat ainds = getainds();
    indexmat dinds = getdinds();

    if (dinds.count() == 0)
        return b.copy();

    vec ba = b.extract(ainds);
    vec bd = b.extract(dinds);

    Vec bapetsc = ba.getpetsc();
    Vec bdpetsc = bd.getpetsc();

    Vec prodvec;
    VecCreate(PETSC_COMM_SELF, &prodvec);
    VecSetSizes(prodvec, PETSC_DECIDE, ainds.count());
    VecSetFromOptions(prodvec);   
    MatMult(getdpetsc(), bdpetsc, prodvec);
    VecAXPY(bapetsc, -1, prodvec);
    VecDestroy(&prodvec);

    return ba;
}
        
indexmat mat::getainds(void) { errorifpointerisnull(); errorifinvalidated(); return rawmatptr->getainds(); }
indexmat mat::getdinds(void) { errorifpointerisnull(); errorifinvalidated(); return rawmatptr->getdinds(); }

Mat mat::getapetsc(void) { errorifpointerisnull(); errorifinvalidated(); return rawmatptr->getapetsc(); }
Mat mat::getdpetsc(void) { errorifpointerisnull(); errorifinvalidated(); return rawmatptr->getdpetsc(); }

void mat::print(void) { errorifpointerisnull(); errorifinvalidated(); rawmatptr->print(); }

mat mat::copy(void)
{
    errorifpointerisnull(); errorifinvalidated();
    
    Mat outa, outd;
    MatDuplicate(getapetsc(), MAT_COPY_VALUES, &outa);
    MatDuplicate(getdpetsc(), MAT_COPY_VALUES, &outd);
    return mat(std::shared_ptr<rawmat>(new rawmat(  rawmatptr->getdofmanager(), outa, outd, getainds().copy(), getdinds().copy()  )));
}


mat mat::operator+(void) { return copy(); }
mat mat::operator-(void) { return *this*-1; }

mat mat::operator*(double input)
{
    errorifpointerisnull(); errorifinvalidated();
    
    Mat outa, outd;
    MatDuplicate(getapetsc(), MAT_DO_NOT_COPY_VALUES, &outa);
    MatDuplicate(getdpetsc(), MAT_DO_NOT_COPY_VALUES, &outd);
    MatAXPY(outa, input, getapetsc(), SAME_NONZERO_PATTERN);
    MatAXPY(outd, input, getdpetsc(), SAME_NONZERO_PATTERN);
    return mat(std::shared_ptr<rawmat>(new rawmat(  rawmatptr->getdofmanager(), outa, outd, getainds().copy(), getdinds().copy()  )));
}

mat mat::operator/(double input) { return *this*(1.0/input); }

mat mat::operator*(mat input)
{
    errorifpointerisnull(); errorifinvalidated();

    // | A   D |  | B   E |   | AB  AE+D |
    // |       |  |       | = |          |
    // | 0   1 |  | 0   1 |   | 0      1 |
    
    Mat A = getapetsc();
    Mat D = getdpetsc();
    Mat B = input.getapetsc();
    Mat E = input.getdpetsc();
    
    Mat AB, AEplusD;
    MatMatMult(A, B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AB);
    MatMatMult(A, E, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AEplusD);
    MatAXPY(AEplusD, 1, D, DIFFERENT_NONZERO_PATTERN);
    
    return mat(std::shared_ptr<rawmat>(new rawmat(  rawmatptr->getdofmanager(), AB, AEplusD, getainds().copy(), getdinds().copy()  )));
}

mat mat::operator+(mat input)
{
    errorifpointerisnull(); errorifinvalidated();
    
    mat copied = copy();
    MatAXPY(copied.getapetsc(), 1, input.getapetsc(), DIFFERENT_NONZERO_PATTERN);
    MatAXPY(copied.getdpetsc(), 1, input.getdpetsc(), DIFFERENT_NONZERO_PATTERN);
    return copied;
}

mat mat::operator-(mat input)
{
    errorifpointerisnull(); errorifinvalidated();
    
    mat copied = copy();
    MatAXPY(copied.getapetsc(), -1, input.getapetsc(), DIFFERENT_NONZERO_PATTERN);
    MatAXPY(copied.getdpetsc(), -1, input.getdpetsc(), DIFFERENT_NONZERO_PATTERN);
    return copied;
}

vec mat::operator*(vec input)
{
    errorifpointerisnull(); errorifinvalidated();
    
    indexmat ainds = getainds();
    indexmat dinds = getdinds();

    vec ia = input.extract(ainds);
    vec id = input.extract(dinds);

    Vec iapetsc = ia.getpetsc();
    Vec idpetsc = id.getpetsc();

    Vec prodvec;
    VecCreate(PETSC_COMM_SELF, &prodvec);
    VecSetSizes(prodvec, PETSC_DECIDE, ainds.count());
    VecSetFromOptions(prodvec);
    MatMult(getapetsc(), iapetsc, prodvec);
    MatMultAdd(getdpetsc(), idpetsc, prodvec, prodvec);

    vec vp(std::shared_ptr<rawvec>(new rawvec(  ia.getpointer()->getdofmanager(), prodvec  )));

    return xbmerge(vp, input);
}


mat operator*(double inputdouble, mat inputmat) { return inputmat*inputdouble; }

