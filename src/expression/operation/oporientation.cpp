#include "oporientation.h"


std::vector<std::vector<densemat>> oporientation::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }

    int numelems = elemselect.countinselection();
    int numevalpts = evaluationcoordinates.size()/3;
    
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    std::vector<bool> prdef = universe::getrawmesh()->getphysicalregions()->get(myphysreg)->getdefinition();
    int elemdim = elemselect.getelementdimension();
    int elemtypenum = elemselect.getelementtypenumber();
    std::vector<int> elemnums = elemselect.getelementnumbers();
    
    std::vector<std::string> geonames = {"point","line","face","volume"};


    if (elemdim != problemdimension-1)
    {
        logs log;
        log.msg() << "Error in 'oporientation' object: operation can only be evaluated on a region of one lower dimension than the geometry" << std::endl;
        log.msg() << "Did you try to compute normal(physreg) on a " << geonames[elemdim] << " region for a " << problemdimension << "D geometry?" << std::endl;
        log.error();
    }
    
    densemat output(numelems,1,1);
    double* outptr = output.getvalues();
    

    // Calculate the orientations:
    std::vector<double>* nodecoords = universe::getrawmesh()->getnodes()->getcoordinates();
    elements* els = universe::getrawmesh()->getelements();
    
    std::vector<int> celltypes(numelems), cellnums(numelems);
    
    bool isnotfound = false, ismultifound = false;
    for (int i = 0; i < numelems; i++)
    {
        // Boundary element number:
        int e = elemnums[i];
        
        // Get the corresponding cell number and type in 'myphysreg' or give an error if not touching:
        std::vector<int> cellstouching = els->getcellsontype(elemtypenum, e);
        
        // Select the cell which is in 'myphysreg':
        int numfound = 0;
        for (int c = 0; c < cellstouching.size()/2; c++)
        {
            int ctt = cellstouching[2*c+0];
            int ctn = cellstouching[2*c+1];
            int curdisjreg = els->getdisjointregion(ctt,ctn);
                
            if (prdef[curdisjreg])
            {
                celltypes[i] = ctt;
                cellnums[i] = ctn;
                numfound++;
            }
        }
        isnotfound = isnotfound || (numfound == 0);
        ismultifound = ismultifound || (numfound > 1);
    }
    
    if (isnotfound)
    {
        logs log;
        log.msg() << "Error in 'oporientation' object: found a " << geonames[elemdim] << " not in contact with the requested physical region" << std::endl;
        log.msg() << "Did you try to compute normal(physreg) on a " << geonames[elemdim] << " not touching the argument region?" << std::endl;
        log.error();
    }
    if (ismultifound)
    {
        logs log;
        log.msg() << "Error in 'oporientation' object: found a " << geonames[elemdim] << " surrounded on both sides by a " << geonames[elemdim+1] << std::endl;
        log.msg() << "Did you try to compute normal(physreg) on a " << geonames[elemdim] << " surrounded on both sides by the argument region?" << std::endl;
        log.error();
    }
    
    indexmat celltypmat(celltypes.size(), 1, celltypes);
    std::vector<std::vector<int>> indexes = celltypmat.findalloccurences(8);
    
    for (int t = 0; t < indexes.size(); t++)
    {
        if (indexes[t].size() == 0)
            continue;
            
        std::vector<int> cellnumsintype;
        gentools::select(cellnums, indexes[t], cellnumsintype);
        
        std::vector<int> subelemnums;
        gentools::select(elemnums, indexes[t], subelemnums);
        
        // 1D must be treated differently:
        if (problemdimension == 1)
        {
            for (int i = 0; i < indexes[t].size(); i++)
            {
                int cornera, cornerb;
                cornera = els->getsubelement(0, 1, cellnumsintype[i], 0);
                if (subelemnums[i] == cornera)
                    cornerb = els->getsubelement(0, 1, cellnumsintype[i], 1);
                else
                {
                    cornerb = cornera;
                    cornera = els->getsubelement(0, 1, cellnumsintype[i], 1);
                }

                if (nodecoords->at(3*cornerb+0) - nodecoords->at(3*cornera+0) > 0)
                    outptr[indexes[t][i]] = -1;
            }
            continue;
        }
        
        // Check if node order flipped compared to in parent:
        std::vector<bool> isflipped = els->isflipped(elemtypenum, subelemnums, t, cellnumsintype);
        
        // Check detjac sign of parent:
        std::vector<int> disjregs = universe::getrawmesh()->getphysicalregions()->get(myphysreg)->getdisjointregionsoftype(t);
        elementselector subselect(disjregs, cellnumsintype, false);
        
        gausspoints gp(t, 0);
        jacobian subjac(subselect, gp.getcoordinates(), NULL);
        densemat subdetjac = subjac.getdetjac();
        
        double* subdetjacptr = subdetjac.getvalues();
        
        for (int i = 0; i < indexes[t].size(); i++)
        {
            bool isflip = isflipped[i];
            bool ispos = (subdetjacptr[i] > 0);
            if (problemdimension == 2)
                ispos = not(ispos);
            
            if (ispos && isflip || not(ispos) && not(isflip))
                outptr[indexes[t][i]] = -1;
        }
    }
    
    output = output.duplicatehorizontally(numevalpts);
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), {{},{output}});
    
    return {{},{output}};
}

densemat oporientation::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densemat output = interpolate(elemselect, evaluationcoordinates, meshdeform)[1][0];
    
    output = output.getflattened();
    output = output.duplicatevertically(numtimeevals);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> oporientation::copy(void)
{
    std::shared_ptr<oporientation> op(new oporientation(myphysreg));
    *op = *this;
    op->reuse = false;
    return op;
}

void oporientation::print(void) {}

