#include "eigenvalue.h"


eigenvalue::eigenvalue(mat A)
{
    myA = A;
    
    // Make sure the Dirichlet constraints have been removed:
    if (A.getpointer()->getdofmanager()->countconstraineddofs() > 0)
    {
        std::cout << "Error in 'eigenvalue' object: remove the Dirichlet constraints in the matrix with .removeconstraints()" << std::endl;
        abort();
    }
}

eigenvalue::eigenvalue(mat A, mat B)
{
    myA = A;
    myB = B;
    
    // Make sure the Dirichlet constraints have been removed:
    if (A.getpointer()->getdofmanager()->countconstraineddofs() > 0 || B.getpointer()->getdofmanager()->countconstraineddofs() > 0)
    {
        std::cout << "Error in 'eigenvalue' object: remove the Dirichlet constraints in the matrix with .removeconstraints()" << std::endl;
        abort();
    }
}

void eigenvalue::compute(int numeigenvaluestocompute, double targeteigenvaluemagnitude)
{
    // Define the slepc eigensolver context:
    EPS eps;
    
    EPSCreate( PETSC_COMM_WORLD, &eps );
    
    // To be general we assume a non-hermitian problem:
    if (myB.getpointer() == NULL)
    {
        EPSSetOperators( eps, myA.getpetsc(), NULL );
        EPSSetProblemType(eps, EPS_NHEP);    
    }
    else
    {
        EPSSetOperators( eps, myA.getpetsc(), myB.getpetsc() );
        EPSSetProblemType(eps, EPS_GNHEP);
    }
    
    // Tell slepc how many eigs we want:
    EPSSetDimensions(eps, numeigenvaluestocompute, PETSC_DECIDE, PETSC_DECIDE);
    
    // Set tolerance and max num of iterations allowed:
    EPSSetTolerances(eps, 1e-6, 100);
    
    // Set the eigenvalue solver:
    EPSSetType(eps, EPSKRYLOVSCHUR);
    
    EPSSetFromOptions(eps);
    
    
    // We use a shift and invert transform:
    ST st;
    EPSGetST(eps, &st);
    STSetType(st, STSINVERT);
    
    
    // We target the eigenvalues with a given magnitude:
     EPSSetTarget(eps, targeteigenvaluemagnitude);
    EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
    
    
    // MUMPS petsc solver context:
    KSP ksp;
    STGetKSP(st, &ksp);
    KSPSetType(ksp, "preonly");
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);
    PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
    
    
    // DO THE ACTUAL RESOLUTION:
    EPSSolve( eps );
    
    // Get the number of eigs found:
    int numeigsfound;
    EPSGetConverged( eps, &numeigsfound );
    
        
    // Get all eigs:
    eigenvaluereal.resize(numeigsfound);
    eigenvalueimaginary.resize(numeigsfound);
    eigenvectorreal.resize(numeigsfound);
    eigenvectorimaginary.resize(numeigsfound);
    
    for (int i = 0; i < numeigsfound; i++)
    {
        double eigvalr, eigvali;
        
        // Create the 'eigvecr' and 'eigveci' vectors based on the dofmanager in myA::
        shared_ptr<rawvec> rawr(new rawvec(myA.getpointer()->getdofmanager())); shared_ptr<rawvec> rawi( new rawvec(myA.getpointer()->getdofmanager()));
        vec eigvecr(rawr); vec eigveci(rawi);
        
        EPSGetEigenpair( eps, i, &eigvalr, &eigvali, eigvecr.getpetsc(), eigveci.getpetsc() );
                
        eigenvaluereal[i] = eigvalr;
        eigenvalueimaginary[i] = eigvali;
        
        eigenvectorreal[i] = eigvecr;
        eigenvectorimaginary[i] = eigveci;
    }
}

std::vector<double> eigenvalue::geteigenvaluerealpart(void) { return eigenvaluereal; }
std::vector<double> eigenvalue::geteigenvalueimaginarypart(void) { return eigenvalueimaginary; }
std::vector<vec> eigenvalue::geteigenvectorrealpart(void) { return eigenvectorreal; }
std::vector<vec> eigenvalue::geteigenvectorimaginarypart(void) { return eigenvectorimaginary; }

void eigenvalue::printeigenvalues(void)
{
    std::cout << std::endl << "Printing the " << count() << " eigenvalues in format real part imaginary part:" << std::endl << std::endl;
    
    for (int i = 0; i < count(); i++)
        std::cout << "#" << std::left << std::setw(7) << i << std::left << std::setw(18) << eigenvaluereal[i] << std::left << std::setw(18) << eigenvalueimaginary[i] << std::endl;
    std::cout << std::endl;
}

void eigenvalue::printeigenfrequencies(void)
{
    std::cout << std::endl << "Printing the " << count() << " eigenfrequencies:" << std::endl << std::endl;
    
    for (int i = 0; i < count(); i++)
        std::cout << "#" << std::left << std::setw(7) << i << std::left << std::setw(18) << std::sqrt(eigenvaluereal[i])/(2*3.14159265359) << " Hz" << std::endl;
    std::cout << std::endl;
}



