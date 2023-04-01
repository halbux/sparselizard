#include "spline.h"
#include "sl.h"
#include "mat.h"
#include "vec.h"
#include "gentools.h"


spline::spline(std::string filename, char delimiter)
{
    std::vector<double> data = sl::loadvector(filename, delimiter, false);
    
    if (data.size()%2 != 0)
    {
        logs log;
        log.msg() << "Error in 'spline' object: expected a vector length multiple of 2 in '" << filename << "' (format {x1,y1,x2,y2,...})" << std::endl;
        log.error();
    }
    
    std::vector<double> xin(data.size()/2);
    std::vector<double> yin(data.size()/2);
    
    for (int i = 0; i < data.size()/2; i++)
    {
        xin[i] = data[2*i+0];
        yin[i] = data[2*i+1];
    }
    
    set(xin,yin);
}

spline::spline(std::vector<double> xin, std::vector<double> yin)
{
    set(xin,yin);
}

void spline::set(std::vector<double>& xin, std::vector<double>& yin)
{
    derivativeorder = 0;
    
    if (xin.size() != yin.size())
    {
        logs log;
        log.msg() << "Error in 'spline' object: x and y dataset sizes do not match" << std::endl;
        log.error();
    }   
    int len = xin.size();
    if (len < 2)
    {
        logs log;
        log.msg() << "Error in 'spline' object: expected at least two data points" << std::endl;
        log.error();
    }   
    
    myx = densemat(len,1);
    myy = densemat(len,1);
    double* xvals = myx.getvalues();
    double* yvals = myy.getvalues();
    
    // First sort ascendingly according to x:
    std::vector<int> reorderingvector;
    gentools::stablesort(0, xin, reorderingvector);
    for (int i = 0; i < len; i++)
    {
        xvals[i] = xin[reorderingvector[i]];
        yvals[i] = yin[reorderingvector[i]];
    }
    xmin = xvals[0]; xmax = xvals[len-1];
    
    double absnoise = noisethreshold*std::abs(xmax-xmin);
    for (int i = 1; i < len; i++)
    {
        if (xvals[i]-xvals[i-1] < absnoise)
        {
            logs log;
            log.msg() << "Error in 'spline' object: distance between two samples is " << (xvals[i]-xvals[i-1]) << " (below noise level " << absnoise << ")" << std::endl;
            log.error();
        }
    }
    
    
    // Create the A matrix and b rhs:
    indexmat Arows(3*len-2,1), Acols(3*len-2,1);
    int* Arowvals = Arows.getvalues();
    int* Acolvals = Acols.getvalues();
    densemat Av(3*len-2,1);
    double* Avals = Av.getvalues();

    densemat bv(len,1);
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
    vec b(len, indexmat(len,1,0,1), bv);
    
    vec k = sl::solve(A,b);
    densemat kv = k.getvalues(indexmat(len,1,0,1));
    double* kvals = kv.getvalues();
    
    
    // Compute the spline parameters a and b:
    mya = densemat(len,1);
    myb = densemat(len,1);
    double* aparamvals = mya.getvalues();
    double* bparamvals = myb.getvalues();
    
    for (int i = 1; i < len; i++)
    {
        aparamvals[i] = kvals[i-1]*(xvals[i]-xvals[i-1])-(yvals[i]-yvals[i-1]);
        bparamvals[i] = -kvals[i]*(xvals[i]-xvals[i-1])+(yvals[i]-yvals[i-1]);
    }
}

spline spline::getderivative(void)
{
    spline spl;
    spl = *this;
    spl.derivativeorder++;

    return spl;
}

double spline::evalat(double input)
{
    std::vector<double> invec = {input};
    return evalat(invec)[0];
}

std::vector<double> spline::evalat(std::vector<double> input)
{
    densemat indm(input.size(),1, input);
    densemat outdm = evalat(indm);
    std::vector<double> output;
    outdm.getvalues(output);
    
    return output;
}

densemat spline::evalat(densemat input)
{
    int numin = input.count();
    double* inputvals = input.getvalues();

    std::vector<double> invals;
    input.getvalues(invals);
    
    // Sort the input data ascendingly:
    std::vector<int> reorderingvector;
    gentools::stablesort(0, invals, reorderingvector);
    for (int i = 0; i < numin; i++)
        inputvals[i] = invals[reorderingvector[i]];
    double inmin = inputvals[0]; double inmax = inputvals[numin-1];
    
    std::vector<double> outvec(numin);
    
    // Get the corresponding data via spline interpolation.
    double* xvals = myx.getvalues(); double* yvals = myy.getvalues();
    double* avals = mya.getvalues(); double* bvals = myb.getvalues();
    
    // First find the corresponding spline:
    int curspline = 1;
    for (int i = 0; i < numin; i++)
    {
        double cur = inputvals[i];
        // Find the spline:
        while (curspline < myx.count()-1 && xvals[curspline] < cur)
            curspline++;
        // Interpolate on the spline:
        double dx = xvals[curspline]-xvals[curspline-1];
        double tx = (cur-xvals[curspline-1])/dx;
        tx = std::max(0.0, std::min(1.0, tx));
        double a = avals[curspline];
        double b = bvals[curspline];
        
        double y = (1.0-tx)*yvals[curspline-1] + tx*yvals[curspline] + tx*(1.0-tx)*((1.0-tx)*a+tx*b);
        double dydx = 1.0/dx * (yvals[curspline]-yvals[curspline-1] + (1.0-2.0*tx)*(a*(1.0-tx)+b*tx) + tx*(1.0-tx)*(b-a));
        
        // The spline turns into a straight line at xmin and xmax:
        if (cur < xmin)
            y = dydx*cur + y - dydx*xmin;
        if (cur > xmax)
            y = dydx*cur + y - dydx*xmax;
        
        if (derivativeorder == 0)
            outvec[i] = y;
        if (derivativeorder == 1)
            outvec[i] = dydx;
        if (derivativeorder == 2)
            outvec[i] = 2.0/(dx*dx) * (b-2.0*a + (a-b)*3.0*tx);
        if (derivativeorder == 3)
            outvec[i] = 6.0/(dx*dx*dx) * (a-b);
        if (derivativeorder > 3)
            outvec[i] = 0.0;
    }
    
    // Unsort the data:
    densemat output(input.countrows(),input.countcolumns());
    double* outputvals = output.getvalues();
    
    for (int i = 0; i < numin; i++)
        outputvals[reorderingvector[i]] = outvec[i];
    
    return output;
}

void spline::write(std::string filename, int numsplits, char delimiter)
{
    if (numsplits < 0)
    {
        logs log;
        log.msg() << "Error in 'spline' object: cannot write with " << numsplits << " splits" << std::endl;
        log.error();
    }

    // Get the x positions:
    double* xvalues = myx.getvalues();
    
    densemat xsplit(1+(myx.count()-1)*(numsplits+1),1);
    double* xsplitvals = xsplit.getvalues();
    double step = 1.0/(numsplits+1.0);
    
    xsplitvals[0] = xvalues[0];
    
    int index = 1;
    for (int i = 0; i < myx.count()-1; i++)
    {
        for (int j = 0; j < numsplits+1; j++)
        {
            xsplitvals[index] = xvalues[i]+(j+1.0)*step*(xvalues[i+1]-xvalues[i]);
            index++;
        }
    }
    densemat evaled = evalat(xsplit);
    double* evaledvals = evaled.getvalues();
    
    
    // Write to file:
    std::vector<double> data(2*xsplit.count());
    for (int i = 0; i < xsplit.count(); i++)
    {
        data[2*i+0] = xsplitvals[i];
        data[2*i+1] = evaledvals[i];
    }
    
    sl::writevector(filename, data, delimiter, false);
}

