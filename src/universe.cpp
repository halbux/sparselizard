#include "universe.h"


double universe::roundoffnoiselevel = 1e-12;

std::shared_ptr<rawmesh> universe::mymesh = NULL;

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

int universe::physregshift = 0;

void universe::allowestimatorupdate(bool allowitonce)
{
    if (allowitonce)
    {
        numallowedtimes++;
        if (numallowedtimes == 1)
            estimatorcalcstate++;
    }
    else
        numallowedtimes--;
}

bool universe::isestimatorupdateallowed(long long int statenumber)
{
    return (numallowedtimes > 0 && statenumber != estimatorcalcstate);
}

long long int universe::estimatorcalcstate = 0;

int universe::numallowedtimes = 0;
        
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

std::tuple<std::shared_ptr<jacobian>, std::vector<std::shared_ptr<operation>>,std::vector<std::shared_ptr<operation>>, std::vector< std::vector<std::vector<densematrix>> >,std::vector< densematrix >> universe::selectsubset(int numevalpts, std::vector<int>& selectedelementindexes)
{
    int numselected = selectedelementindexes.size();

    auto output = std::make_tuple(computedjacobian, oppointers, oppointersfft, opcomputed, opcomputedfft);
    
    // Replace the jacobian with a subset of it:
    if (computedjacobian != NULL)
    {
        std::shared_ptr<jacobian> newjac(new jacobian);
        *newjac = computedjacobian->extractsubset(selectedelementindexes);
        computedjacobian = newjac;
    }
        
    // Replace the computed ops by their subset:
    for (size_t i = 0; i < opcomputed.size(); i++)
    {
        for (size_t h = 0; h < opcomputed[i].size(); h++)
        {
            if (opcomputed[i][h].size() > 0)
                opcomputed[i][h][0] = opcomputed[i][h][0].extractrows(selectedelementindexes);
        }
    }
    
    // In the multiharmonic case columns must be selected:
    std::vector<int> mhcols;
    if (opcomputedfft.size() > 0)
    {
        mhcols = std::vector<int>(numselected*numevalpts);
        for (int i = 0; i < numselected; i++)
        {
            for (int j = 0; j < numevalpts; j++)
                mhcols[i*numevalpts+j] = selectedelementindexes[i]*numevalpts+j;
        }
    }
    for (size_t i = 0; i < opcomputedfft.size(); i++)
        opcomputedfft[i] = opcomputedfft[i].extractcols(mhcols);
    
    return output;
}

void universe::restore(std::tuple<std::shared_ptr<jacobian>, std::vector<std::shared_ptr<operation>>,std::vector<std::shared_ptr<operation>>, std::vector< std::vector<std::vector<densematrix>> >,std::vector< densematrix >> input)
{
    computedjacobian = std::get<0>(input);
    
    oppointers = std::get<1>(input);
    oppointersfft = std::get<2>(input);
    opcomputed = std::get<3>(input);
    opcomputedfft = std::get<4>(input);
}


std::shared_ptr<jacobian> universe::computedjacobian = NULL;


std::vector<std::shared_ptr<operation>> universe::oppointers = {};
std::vector<std::shared_ptr<operation>> universe::oppointersfft = {};
std::vector< std::vector<std::vector<densematrix>> > universe::opcomputed = {};
std::vector< densematrix > universe::opcomputedfft = {};

int universe::getindexofprecomputedvalue(std::shared_ptr<operation> op)
{
    for (size_t i = 0; i < oppointers.size(); i++)
    {
        if (oppointers[i].get() == op.get())
            return i;
    }
    return -1;
}

int universe::getindexofprecomputedvaluefft(std::shared_ptr<operation> op)
{
    for (size_t i = 0; i < oppointersfft.size(); i++)
    {
        if (oppointersfft[i].get() == op.get())
            return i;
    }
    return -1;
}

int universe::getindexofprecomputedvalue(std::shared_ptr<rawparameter> param, int row, int col)
{
    for (size_t i = 0; i < oppointers.size(); i++)
    {
        if (oppointers[i]->isparameter() && (oppointers[i]->getparameterpointer()).get() == param.get() && oppointers[i]->getselectedrow() == row && oppointers[i]->getselectedcol() == col)
            return i;
    }
    return -1;
}

int universe::getindexofprecomputedvaluefft(std::shared_ptr<rawparameter> param, int row, int col)
{
    for (size_t i = 0; i < oppointersfft.size(); i++)
    {
        if (oppointersfft[i]->isparameter() && (oppointersfft[i]->getparameterpointer()).get() == param.get() && oppointersfft[i]->getselectedrow() == row && oppointersfft[i]->getselectedcol() == col)
            return i;
    }
    return -1;
}

int universe::getindexofprecomputedvalue(std::shared_ptr<rawfield> rf, int td, int sd, int kepd, int ffc)
{
    for (size_t i = 0; i < oppointers.size(); i++)
    {
        if (oppointers[i]->isfield() && (oppointers[i]->getfieldpointer()).get() == rf.get() && oppointers[i]->getformfunctioncomponent() == ffc && oppointers[i]->getspacederivative() == sd && oppointers[i]->getkietaphiderivative() == kepd && oppointers[i]->gettimederivative() == td)
            return i;
    }
    return -1;
}

int universe::getindexofprecomputedvaluefft(std::shared_ptr<rawfield> rf, int td, int sd, int kepd, int ffc)
{
    for (size_t i = 0; i < oppointersfft.size(); i++)
    {
        if (oppointersfft[i]->isfield() && (oppointersfft[i]->getfieldpointer()).get() == rf.get() && oppointersfft[i]->getformfunctioncomponent() == ffc && oppointersfft[i]->getspacederivative() == sd && oppointersfft[i]->getkietaphiderivative() == kepd && oppointersfft[i]->gettimederivative() == td)
            return i;
    }
    return -1;
}

std::vector<std::vector<densematrix>> universe::getprecomputed(int index)
{
    std::vector<std::vector<densematrix>> output = opcomputed[index];
    for (size_t h = 0; h < output.size(); h++)
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

void universe::setprecomputed(std::shared_ptr<operation> op, std::vector<std::vector<densematrix>> val)
{
    oppointers.push_back(op);
    opcomputed.push_back(val);
    for (size_t h = 0; h < val.size(); h++)
    {
        if (val[h].size() == 1)
            opcomputed[opcomputed.size()-1][h][0] = val[h][0].copy();
    }
}

void universe::setprecomputedfft(std::shared_ptr<operation> op, densematrix val)
{
    oppointersfft.push_back(op);
    opcomputedfft.push_back(val.copy());
}

bool universe::keeptrackofrhsassembly = false;
std::vector<std::pair<intdensematrix, densematrix>> universe::rhsterms = {}; 
        
std::vector<std::vector<vec>> universe::xdtxdtdtx = {{},{},{}};        



std::vector<std::pair< std::string, std::vector<std::vector< std::vector<hierarchicalformfunctioncontainer> >> >> universe::formfuncpolys = {};     

hierarchicalformfunctioncontainer* universe::gethff(std::string fftypename, int elementtypenumber, int interpolorder, std::vector<double> evaluationcoordinates)
{
    // Find the type name in the container:
    int typenameindex = -1;
    for (size_t i = 0; i < formfuncpolys.size(); i++)
    {
        if (formfuncpolys[i].first == fftypename)
        {
            typenameindex = i;
            break;
        }
    }

    // In case the form function polynomials are available:
    if (typenameindex != -1 && formfuncpolys[typenameindex].second[elementtypenumber].size() > (size_t) interpolorder &&
        formfuncpolys[typenameindex].second[elementtypenumber][interpolorder].size() > 0)
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
    if (formfuncpolys[typenameindex].second[elementtypenumber].size() <= (size_t) interpolorder)
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
    for (size_t typenameindex = 0; typenameindex < formfuncpolys.size(); typenameindex++)
    {
        for (size_t elementtypenumber = 0; elementtypenumber < (formfuncpolys[typenameindex].second).size(); elementtypenumber++)
        {
            for (size_t interpolorder = 0; interpolorder < formfuncpolys[typenameindex].second[elementtypenumber].size(); interpolorder++)
            {
                if (formfuncpolys[typenameindex].second[elementtypenumber][interpolorder].size() > 0)
                    formfuncpolys[typenameindex].second[elementtypenumber][interpolorder][0].setvaluestatus(false);
            }
        }
    }
}


std::vector<std::vector<std::vector<std::vector<int>>>> universe::splitdefinition = std::vector<std::vector<std::vector<std::vector<int>>>>(8, std::vector<std::vector<std::vector<int>>>(0));

bool universe::getsplitdefinition(std::vector<std::vector<int>>& splitdef, int elementtypenumber, int splitnum, std::vector<int>& edgenumbers)
{
    int ne = edgenumbers.size();
    int numrel = myalgorithm::factorial(ne);
    int rel = myalgorithm::identifyrelations(edgenumbers);
    
    if (splitdefinition[elementtypenumber].size() == 0 || splitdefinition[elementtypenumber][splitnum*numrel+rel].size() == 0)
        return false;
        
    splitdef = splitdefinition[elementtypenumber][splitnum*numrel+rel];
    return true;
}

void universe::setsplitdefinition(std::vector<std::vector<int>>& splitdef, int elementtypenumber, int splitnum, std::vector<int>& edgenumbers)
{
    int ne = edgenumbers.size();
    int numrel = myalgorithm::factorial(ne);
    int rel = myalgorithm::identifyrelations(edgenumbers);
    
    if (splitdefinition[elementtypenumber].size() == 0)
        splitdefinition[elementtypenumber] = std::vector<std::vector<std::vector<int>>>(std::pow(2,ne)*numrel, std::vector<std::vector<int>>(0));
    
    splitdefinition[elementtypenumber][splitnum*numrel+rel] = splitdef;
}

