#include "petscmesh.h"
#include "universe.h"


void petscmesh::reordernodes(int ourtypenum, std::vector<int>& toreorder)
{
    switch (ourtypenum)
    {
        case 5:
        {
            toreorder = {toreorder[0], toreorder[1], toreorder[2], toreorder[3], toreorder[4], toreorder[7], toreorder[6], toreorder[5]};
            break;
        }
        case 6:
        {
            toreorder = {toreorder[0], toreorder[2], toreorder[1], toreorder[3], toreorder[4], toreorder[5]};
            break;
        }
    }
}

petscmesh::petscmesh(std::string filename)
{
    std::vector<std::string> supportedextensions = {".msh",".cgns",".exo",".gen",".cas",".h5",".med",".ply",".dat"};
 
    bool isvalidext = false;
    
    std::string curfileext = myalgorithm::getfileextension(filename);
    for (int i = 0; i < supportedextensions.size(); i++)
    {
        if (supportedextensions[i] == curfileext)
        {
            isvalidext = true;
            break;
        }
    }
    
    if (isvalidext)
    {
        DMPlexCreateFromFile(PETSC_COMM_SELF, filename.c_str(), PETSC_TRUE, &mypetscmesh);
        DMGetDimension(mypetscmesh, &meshdim);
        return;
    }
    
    std::cout << "Error in 'petscmesh' object: cannot load mesh file '" << filename << "'." << std::endl;
    std::cout << "Supported mesh formats are ";
    for (int i = 0; i < supportedextensions.size()-1; i++)
        std::cout << "'" << supportedextensions[i] << "', ";
    std::cout << "'" << supportedextensions[supportedextensions.size()-1] << "'." << std::endl;
    abort();
}

petscmesh::~petscmesh(void)
{
    DMDestroy(&mypetscmesh);
}

void petscmesh::extract(nodes& mynodes, elements& myelements, physicalregions& myphysicalregions, bool verbosity)
{
    ///// Extract the node coordinates:
    
    Vec coordvec = PETSC_NULL;
    DMGetCoordinates(mypetscmesh, &coordvec);
    
    int numberofnodes;
    VecGetSize(coordvec, &numberofnodes);
    numberofnodes = numberofnodes/meshdim;
    
    mynodes.setnumber(numberofnodes);
    // Transfer from Vec to the node object:
    intdensematrix addresses(meshdim*numberofnodes,1, 0,1);
    densematrix coordmat(meshdim*numberofnodes,1);
    VecGetValues(coordvec, meshdim*numberofnodes, addresses.getvalues(), coordmat.getvalues());
    
    double* coordmatval = coordmat.getvalues();
    std::vector<double>* nodecoordinates = mynodes.getcoordinates();
    for (int i = 0; i < numberofnodes; i++)
    {
        for (int j = 0; j < meshdim; j++)
            nodecoordinates->at(3*i+j) = coordmatval[meshdim*i+j];
    }
    
    if (verbosity)
        mynodes.print();
    
    
    ///// Extract the elements and the physical regions:
    
    // Point range for every dimension:
    std::vector<int> rangestart(meshdim+1), rangeend(meshdim+1);
    for (int i = 0; i <= meshdim; i++)
        DMPlexGetDepthStratum(mypetscmesh, i, &rangestart[i], &rangeend[i]);
    // Vertex number range [vstart,vend[
    int vstart = rangestart[0], vend = rangeend[0];
    
    // Number of labels:
    int numlabels;
    DMGetNumLabels(mypetscmesh, &numlabels);
    
    element myelem(0);

    // The first two labels are for the cell type and the depth (constructed automatically):
    for (int l = 2; l < numlabels; l++)
    {
        DMLabel curlabel;
        DMGetLabelByNum(mypetscmesh, l, &curlabel);
        
        int numphysregsinlabel;
        DMLabelGetNumValues(curlabel, &numphysregsinlabel);
        
        IS labelIS;
        DMLabelGetValueIS(curlabel, &labelIS);

        const int* physregsinlabel = NULL;
        ISGetIndices(labelIS, &physregsinlabel);

        for (int p = 0; p < numphysregsinlabel; p++)
        {
            int curphysreg = physregsinlabel[p];

            IS pointIS;
            DMLabelGetStratumIS(curlabel, curphysreg, &pointIS);

            int numptsincurphysreg;
            ISGetSize(pointIS, &numptsincurphysreg);
            
            const int* curpoints = NULL;
            ISGetIndices(pointIS, &curpoints);
            
            for (int i = 0; i < numptsincurphysreg; i++)
            {
                int curpt = curpoints[i];
                
                // Get the dimension of the current element:
                int dim = -1;
                for (int j = 0; j <= meshdim; j++)
                {
                    if (curpt >= rangestart[j] && curpt < rangeend[j])
                    {
                        dim = j; 
                        break;
                    }
                }
             
                int curlen;
                int* points = NULL;
                DMPlexGetTransitiveClosure(mypetscmesh, curpt, PETSC_TRUE, &curlen, &points);

                // Count number of nodes in the element:
                int numnodesinelem = 0;
                for (int j = 0; j < 2*curlen; j=j+2)
                {
                    int cur = points[j];
                    if (cur >= vstart && cur < vend)
                        numnodesinelem++;
                }
                // Deduce the element type:
                int elemtypenum = myelem.deducetypenumber(dim, numnodesinelem);
                // Put the nodes in a vector:
                std::vector<int> nodesinelem(numnodesinelem);
                int index = 0;
                for (int j = 0; j < 2*curlen; j=j+2)
                {
                    int cur = points[j];
                    if (cur >= vstart && cur < vend)
                    {
                        nodesinelem[index] = cur - vstart; 
                        index++;
                    }
                }
                
                DMPlexRestoreTransitiveClosure(mypetscmesh, curpt, PETSC_TRUE, &curlen, &points);
             
                // Bring to our element node ordering:
                reordernodes(elemtypenum, nodesinelem);
                
                // Add element:
                int elementindexincurrenttype = myelements.add(elemtypenum, curvatureorder, nodesinelem);
                
                // Get the physical region and add the element:
                physicalregion* currentphysicalregion = myphysicalregions.get(universe::physregshift*(dim+1) + curphysreg);
                currentphysicalregion->addelement(elemtypenum, elementindexincurrenttype);
            }
                
            ISRestoreIndices(pointIS, &curpoints);
            ISDestroy(&pointIS);
        }
        
        ISRestoreIndices(labelIS, &physregsinlabel);
        ISDestroy(&labelIS);
    }
    
    if (verbosity)
    {
        myelements.printnumber();
        myelements.printsubelements();   
    }
}

void petscmesh::view(void)
{
    PetscViewer v;
    PetscViewerFormat format = PETSC_VIEWER_ASCII_CSV;
    PetscViewerASCIIOpen(PETSC_COMM_SELF, "petscview.csv", &v);
    PetscViewerPushFormat(v, format);
    DMView(mypetscmesh, v);
    PetscViewerDestroy(&v);
}

