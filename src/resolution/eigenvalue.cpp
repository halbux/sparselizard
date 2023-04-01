#include "eigenvalue.h"
#include <slepceps.h>
#include <slepcpep.h>


eigenvalue::eigenvalue(mat A)
{
    myA = A;
}

eigenvalue::eigenvalue(mat A, mat B)
{
    myA = A;
    myB = B;
}

eigenvalue::eigenvalue(mat K, mat C, mat M)
{
    mymats = {K,C,M};
}

eigenvalue::eigenvalue(std::vector<mat> inmats)
{
    if (inmats.size() == 0)
    {
        logs log;
        log.msg() << "Error in 'eigenvalue' object: expected at least one matrix in the mat vector" << std::endl;
        log.error();
    }

    mymats = inmats;
}

void eigenvalue::compute(int numeigenvaluestocompute, double targeteigenvaluemagnitude)
{
    if (mymats.size() == 0)
    {
        // Define the slepc eigensolver context:
        EPS eps;
        
        EPSCreate( PETSC_COMM_SELF, &eps );
        
        // To be general we assume a non-hermitian problem:
        if (myB.getpointer() == NULL)
        {
            EPSSetOperators( eps, myA.getapetsc(), NULL );
            EPSSetProblemType(eps, EPS_NHEP);    
        }
        else
        {
            EPSSetOperators( eps, myA.getapetsc(), myB.getapetsc() );
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
        PCFactorSetMatSolverType(pc, universe::solvertype);
        
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
            
            // Create the 'eigvecr' and 'eigveci' vectors:
            std::shared_ptr<rawvec> rawr(new rawvec(std::shared_ptr<dofmanager>(new dofmanager(myA.getainds().count())))); 
            std::shared_ptr<rawvec> rawi(new rawvec(std::shared_ptr<dofmanager>(new dofmanager(myA.getainds().count())))); 
            vec eigvecr(rawr); vec eigveci(rawi);
            
            EPSGetEigenpair( eps, i, &eigvalr, &eigvali, eigvecr.getpetsc(), eigveci.getpetsc() );
                    
            eigenvaluereal[i] = eigvalr;
            eigenvalueimaginary[i] = eigvali;
            
            eigenvectorreal[i] = myA.x0merge(eigvecr);
            eigenvectorimaginary[i] = myA.x0merge(eigveci);
        }
    }
    else
    {
        // Define the slepc eigensolver context:
        PEP pep;
        
        PEPCreate( PETSC_COMM_SELF, &pep );
        
        Mat* petscmats = new Mat[mymats.size()];
        for (int i = 0; i < mymats.size(); i++)
            petscmats[i] = mymats[i].getapetsc();
        
        PEPSetOperators( pep, mymats.size(), petscmats );
        // We assume a general problem:
        PEPSetProblemType( pep, PEP_GENERAL);
        // Tell slepc how many eigs we want:
        PEPSetDimensions( pep, numeigenvaluestocompute, PETSC_DECIDE, PETSC_DECIDE);
        // Set tolerance and max num of iterations allowed:
        PEPSetTolerances( pep, 1e-6, 100);
        // Set the eigenvalue solver:
        PEPSetType( pep, PEPTOAR);

        // We use a shift and invert transform:
        ST st;
        PEPGetST(pep, &st);
        STSetType(st, STSINVERT);
        
        // We target the eigenvalues with a given magnitude:
        PEPSetTarget(pep, targeteigenvaluemagnitude);
        PEPSetWhichEigenpairs(pep, PEP_TARGET_MAGNITUDE);

        // MUMPS petsc solver context:
        KSP ksp;
        STGetKSP(st, &ksp);
        KSPSetType(ksp, "preonly");
        PC pc;
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCLU);

        PEPSTOARSetDetectZeros(pep,PETSC_TRUE);
        PEPSetScale(pep, PEP_SCALE_SCALAR, PETSC_DECIDE, PETSC_NULL, PETSC_NULL, PETSC_DECIDE, PETSC_DECIDE);

        PCFactorSetMatSolverType(pc, universe::solvertype);
        PEPSetFromOptions(pep);
        PEPSetUp(pep);

        // DO THE ACTUAL RESOLUTION:
        PEPSolve( pep );
        // Get the number of eigs found:
        int numeigsfound;
        PEPGetConverged( pep, &numeigsfound );


        // Get all eigs:
        eigenvaluereal.resize(numeigsfound);
        eigenvalueimaginary.resize(numeigsfound);
        eigenvectorreal.resize(numeigsfound);
        eigenvectorimaginary.resize(numeigsfound);
        
        for (int i = 0; i < numeigsfound; i++)
        {
            double eigvalr, eigvali;
            
            // Create the 'eigvecr' and 'eigveci' vectors:
            std::shared_ptr<rawvec> rawr(new rawvec(std::shared_ptr<dofmanager>(new dofmanager(mymats[0].getainds().count())))); 
            std::shared_ptr<rawvec> rawi(new rawvec(std::shared_ptr<dofmanager>(new dofmanager(mymats[0].getainds().count())))); 
            vec eigvecr(rawr); vec eigveci(rawi);
            
            PEPGetEigenpair( pep, i, &eigvalr, &eigvali, eigvecr.getpetsc(), eigveci.getpetsc() );
                    
            eigenvaluereal[i] = eigvalr;
            eigenvalueimaginary[i] = eigvali;
            
            eigenvectorreal[i] = mymats[0].x0merge(eigvecr);
            eigenvectorimaginary[i] = mymats[0].x0merge(eigveci);
        }
        
        delete[] petscmats;
    }
}

int eigenvalue::count(void) { return eigenvaluereal.size(); }

std::vector<double> eigenvalue::geteigenvaluerealpart(void) { return eigenvaluereal; }
std::vector<double> eigenvalue::geteigenvalueimaginarypart(void) { return eigenvalueimaginary; }
std::vector<vec> eigenvalue::geteigenvectorrealpart(void) { return eigenvectorreal; }
std::vector<vec> eigenvalue::geteigenvectorimaginarypart(void) { return eigenvectorimaginary; }

void eigenvalue::printeigenvalues(void)
{
    std::cout << std::endl << "Printing the " << count() << " eigenvalues in format real part imaginary part:" << std::endl << std::endl;
    
    for (int i = 0; i < count(); i++)
        std::cout << "#" << std::left << std::setw(7) << i << std::left << std::setw(16) << eigenvaluereal[i] << std::left << std::setw(16) << eigenvalueimaginary[i] << std::endl;
    std::cout << std::endl;
}

void eigenvalue::printeigenfrequencies(void)
{   
    if (count() == 0)
        return;

    double pival = 3.1415926535897932384;

    if (mymats.size() == 0 && myB.getpointer() != NULL)
    {
        std::cout << std::endl << "Printing the " << count() << " eigenfrequencies [Hz]:" << std::endl << std::endl;
    
        for (int i = 0; i < count(); i++)
        {
            std::cout << "#" << std::left << std::setw(7) << i+1 << std::left << std::setw(16) << std::sqrt(eigenvaluereal[i])/(2.0*pival);
            if (eigenvalueimaginary[i] != 0)
                std::cout << "Expected a real eigenvalue";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    if (mymats.size() == 3)
    {
        std::cout << std::endl << "Printing the " << count() << " eigenfrequencies:" << std::endl << std::endl;

        std::cout << "        ";
        std::cout << std::left << std::setw(17) << "Damped [Hz]";
        std::cout << std::left << std::setw(17) << "Undamped* [Hz]";
        std::cout << std::left << std::setw(17) << "Bandwidth** [Hz]";
        std::cout << std::left << std::setw(17) << "Damping ratio";
        std::cout << std::left << std::setw(17) << "Q factor" << std::endl;

        double fd,fud,zeta,imf,Qf,Bw;

        for (int i = 0; i < count(); i++)
        {	
            std::cout << "#" << std::left << std::setw(7) << i+1;

            fd = std::abs(eigenvalueimaginary[i])/(2.0*pival);
            imf = std::abs(eigenvaluereal[i])/(2.0*pival);
            fud = std::sqrt( std::pow(eigenvaluereal[i],2.0) + std::pow(eigenvalueimaginary[i],2.0) ) ;
            zeta = -eigenvaluereal[i]/fud;
            fud = fud/(2.0*pival);
            Qf = fud/(2.0*imf);
            Bw = fd/Qf;

            std::cout << std::left << std::setw(17) << fd;
            std::cout << std::left << std::setw(17) << fud;
            std::cout << std::left << std::setw(17) << Bw;
            std::cout << std::left << std::setw(17) << zeta;
            std::cout << std::left << std::setw(17) << Qf << std::endl;
        }
        std::cout << std::endl;
        std::cout << "*Only valid for proportional damping" << std::endl;
        std::cout << "**At -3dB (71% of peak signal)" << std::endl;
        std::cout << std::endl;
    }
}

