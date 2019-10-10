#include "splines.h"


splines::splines(std::vector<double>& xin, std::vector<double>& yin)
{
    if (xin.size() != yin.size())
    {
        std::cout << "Error in 'spline' object: x and y dataset sizes do not match" << std::endl;
        abort();
    }   
    int len = xin.size();
    if (len < 2)
    {
        std::cout << "Error in 'spline' object: expected at least two data points" << std::endl;
        abort();
    }   
    
    myx = densematrix(len,1);
    myy = densematrix(len,1);
    double* xvals = myx.getvalues();
    double* yvals = myy.getvalues();
    
    // First sort ascendingly according to x:
    std::vector<int> reorderingvector;
    myalgorithm::stablesort(0, xin, reorderingvector);
    for (int i = 0; i < reorderingvector.size(); i++)
    {
        xvals[i] = xin[reorderingvector[i]];
        yvals[i] = yin[reorderingvector[i]];
    }
    
    // Create the A matrix and b rhs:
    intdensematrix Arows(3*len-2,1), Acols(3*len-2,1);
    int* Arowvals = Arows.getvalues();
    int* Acolvals = Acols.getvalues();
    densematrix Av(3*len-2,1);
    double* Avals = Av.getvalues();

    densematrix bv(len,1);
    double* bvals = bv.getvalues();
    
    // First row has no left neighbour:
    Arowvals[0] = 0; Arowvals[1] = 0; Acolvals[0] = 0; Acolvals[1] = 1;
    Avals[0] = 2.0/(xvals[1]-xvals[0]); Avals[1] = 1.0/(xvals[1]-xvals[0]);
    bvals[0] = 3.0*(yvals[1]-yvals[0])*Avals[1]*Avals[1];
    // Rows with two neighbours:
    double b1 = bvals[0]; int ind = 2;
    for (int i = 1; i < len-1; i++)
    {
        Arowvals[ind+0] = i; Arowvals[ind+1] = i; Arowvals[ind+2] = i;
        Acolvals[ind+0] = i-1; Acolvals[ind+1] = i; Acolvals[ind+2] = i+1;
        
        Avals[ind+0] = Avals[ind-1]; Avals[ind+2] = 1.0/(xvals[i+1]-xvals[i]); Avals[ind+1] = 2.0*(Avals[ind+0]+Avals[ind+2]);
        double b2 = 3.0*(yvals[i+1]-yvals[i])*Avals[ind+2]*Avals[ind+2];
        bvals[i] = b1+b2;
        b1 = b2;
        ind += 3;
    }
    // Last row has no right neighbour:
    Arowvals[ind+0] = len-1; Arowvals[ind+1] = len-1; Acolvals[ind+0] = len-2; Acolvals[ind+1] = len-1;
    Avals[ind+0] = 1.0/(xvals[len-1]-xvals[len-2]); Avals[ind+1] = 2.0/(xvals[len-1]-xvals[len-2]);
    bvals[len-1] = 3.0*(yvals[len-1]-yvals[len-2])*Avals[ind+0]*Avals[ind+0];
    
    
    // Solve problem Ak = b:
    mat A(len, Arows, Acols, Av);
    vec b(len, intdensematrix(len,1,0,1), bv);
    
    vec k = mathop::solve(A,b);
    densematrix kv = k.getvalues(intdensematrix(len,1,0,1));
    
    
    densematrix myy0,myy1,mya,myb;
    
}

densematrix splines::evalat(densematrix input)
{
        //errorifoutofrange = outofrangeerror;//////////////////NO! ALWAYS GIVE OUT OF RANGE ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

}

void splines::write(std::string filename, int numevalsperspline)
{

}

