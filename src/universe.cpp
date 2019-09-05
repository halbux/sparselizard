#include "universe.h"


mesh* universe::mymesh = NULL;

bool universe::isaxisymmetric = false;

double universe::currenttimestep = 0;

double universe::fundamentalfrequency = -1;
double universe::getfundamentalfrequency(void)
{
    if (fundamentalfrequency > 0)
        return fundamentalfrequency;
    else
    {
        std::cout << "Error in 'universe' object: the fundamental frequency cannot be negative or 0 (make sure it was set)" << std::endl;
        abort();
    }
}

bool universe::isreuseallowed = false;

void universe::allowreuse(void)
{
    isreuseallowed = true;
}

void universe::forbidreuse(void)
{
    isreuseallowed = false;
    
    resethff();
    
    computedjacobian = NULL;
    
    oppointers = {};
    oppointersfft = {};
    opcomputed = {};
    opcomputedfft = {};
}


shared_ptr<jacobian> universe::computedjacobian = NULL;


std::vector<shared_ptr<operation>> universe::oppointers = {};
std::vector<shared_ptr<operation>> universe::oppointersfft = {};
std::vector< std::vector<std::vector<densematrix>> > universe::opcomputed = {};
std::vector< densematrix > universe::opcomputedfft = {};

int universe::getindexofprecomputedvalue(shared_ptr<operation> op)
{
    for (int i = 0; i < oppointers.size(); i++)
    {
        if (oppointers[i].get() == op.get())
            return i;
    }
    return -1;
}

int universe::getindexofprecomputedvaluefft(shared_ptr<operation> op)
{
    for (int i = 0; i < oppointersfft.size(); i++)
    {
        if (oppointersfft[i].get() == op.get())
            return i;
    }
    return -1;
}

std::vector<std::vector<densematrix>> universe::getprecomputed(int index)
{
    std::vector<std::vector<densematrix>> output = opcomputed[index];
    for (int h = 0; h < output.size(); h++)
    {
        if (output[h].size() == 1)
            output[h][0] = output[h][0].copy();
    }
    return output;
}

densematrix universe::getprecomputedfft(int index)
{
    return (opcomputedfft[index]).copy();
}

void universe::setprecomputed(shared_ptr<operation> op, std::vector<std::vector<densematrix>> val)
{
    oppointers.push_back(op);
    opcomputed.push_back(val);
    for (int h = 0; h < val.size(); h++)
    {
        if (val[h].size() == 1)
            opcomputed[opcomputed.size()-1][h][0] = val[h][0].copy();
    }
}

void universe::setprecomputedfft(shared_ptr<operation> op, densematrix val)
{
    oppointersfft.push_back(op);
    opcomputedfft.push_back(val.copy());
}
        
bool universe::skipgausspointweightproduct = false;

bool universe::keeptrackofrhsassembly = false;
std::vector<std::pair<intdensematrix, densematrix>> universe::rhsterms = {}; 
        
std::vector<std::vector<vec>> universe::xdtxdtdtx = {{},{},{}};        



std::vector<std::pair< std::string, std::vector<std::vector< std::vector<hierarchicalformfunctioncontainer> >> >> universe::formfuncpolys = {};     

hierarchicalformfunctioncontainer* universe::gethff(std::string fftypename, int elementtypenumber, int interpolorder, std::vector<double> evaluationcoordinates)
{
    // Find the type name in the container:
    int typenameindex = -1;
    for (int i = 0; i < formfuncpolys.size(); i++)
    {
        if (formfuncpolys[i].first == fftypename)
        {
            typenameindex = i;
            break;
        }
    }

    // In case the form function polynomials are available:
    if (typenameindex != -1 && formfuncpolys[typenameindex].second[elementtypenumber].size() > interpolorder && formfuncpolys[typenameindex].second[elementtypenumber][interpolorder].size() > 0)
    {
        if (isreuseallowed && formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0].isvalueready())
            return &(formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0]);
        else
        {
            formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0].evaluate(evaluationcoordinates);
            if (isreuseallowed)
                formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0].setvaluestatus(true);
            return &(formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0]);
        }
    }
    
    // Otherwise compute the form function polynomials and store them:
    if (typenameindex == -1)
    {
        formfuncpolys.push_back( std::make_pair(fftypename, std::vector<std::vector< std::vector<hierarchicalformfunctioncontainer> >>(8,std::vector< std::vector<hierarchicalformfunctioncontainer> >(0))) );
        typenameindex = formfuncpolys.size() - 1;
    }
    if (formfuncpolys[typenameindex].second[elementtypenumber].size() <= interpolorder)
        formfuncpolys[typenameindex].second[elementtypenumber].resize(interpolorder+1);

    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, fftypename);
    
    formfuncpolys[typenameindex].second[elementtypenumber][interpolorder] = {myformfunction->evalat(interpolorder)};
    formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0].evaluate(evaluationcoordinates);
    
    if (isreuseallowed)
        formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0].setvaluestatus(true);
    
    return &(formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0]);
}

void universe::resethff(void)
{
    for (int typenameindex = 0; typenameindex < formfuncpolys.size(); typenameindex++)
    {
        for (int elementtypenumber = 0; elementtypenumber < (formfuncpolys[typenameindex].second).size(); elementtypenumber++)
        {
            for (int interpolorder = 0; interpolorder < formfuncpolys[typenameindex].second[elementtypenumber].size(); interpolorder++)
            {
                if (formfuncpolys[typenameindex].second[elementtypenumber][interpolorder].size() > 0)
                    formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0].setvaluestatus(false);
            }
        }
    }
}

