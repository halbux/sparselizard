#include "petscmesh.h"


void petscmesh::reordernodes(int ourtypenum, std::vector<int>& toreorder)
{
    switch (ourtypenum)
    {
        case 5:
            toreorder = {toreorder[0], toreorder[1], toreorder[2], toreorder[3], toreorder[4], toreorder[7], toreorder[6], toreorder[5]};
            break;
    }
}

petscmesh::petscmesh(std::string filename)
{
    // Make sure the format is supported (and has been validated):
    if (filename.size() >= 5)
    {
        // Get the extension:
        std::string fileext = filename.substr(filename.size()-4,4);

        if (fileext == ".exo" || fileext == ".med" || fileext == ".msh" || fileext == ".ply")
        {
            DMPlexCreateFromFile(PETSC_COMM_SELF, filename.c_str(), PETSC_TRUE, &mypetscmesh);
            DMGetDimension(mypetscmesh, &meshdim);
            return;
        }
    }
    
    std::cout << "Error while loading mesh file '" << filename << "'." << std::endl;
    std::cout << "Supported mesh formats are .exo (ExodusII), .med (Salome), .msh (Fluent), .msh (Gmsh) and .ply (Polygon)." << std::endl;
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
    
    int numpoints = 0;
    int* points = NULL;
    
    // Vertex number range [vstart,vend[
    int vstart, vend;
    DMPlexGetDepthStratum(mypetscmesh, 0, &vstart, &vend);
    
    // Number of labels:
    int numlabels;
    DMGetNumLabels(mypetscmesh, &numlabels);
    // The last label is the depth label (constructed automatically):
    numlabels--;
    
    element myelem(0);
    
    // Loop on all objects (vertices, edges, faces and volumes):
    for (int dim = 0; dim <= meshdim; dim++)
    {
        // Number of objects of current dimension:
        int rangestart, rangeend;
        DMPlexGetDepthStratum(mypetscmesh, dim, &rangestart, &rangeend);
        
        for (int l = 0; l < numlabels; l++)
        {
            DMLabel curlabel;
            DMGetLabelByNum(mypetscmesh, l, &curlabel);
            // Value of the label for the object in the loop below:
            int curlabelvalue;

            for (int i = rangestart; i < rangeend; i++)
            {
                DMLabelGetValue(curlabel, i, &curlabelvalue);

                // Skip elements that have no label:
                if (curlabelvalue == -1)
                    continue;
                
                DMPlexGetTransitiveClosure(mypetscmesh, i, PETSC_TRUE, &numpoints, &points);
                
                // Count number of nodes in the element:
                int numnodesinelem = 0;
                for (int j = 0; j < 2*numpoints; j=j+2)
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
                for (int j = 0; j < 2*numpoints; j=j+2)
                {
                    int cur = points[j];
                    if (cur >= vstart && cur < vend)
                    {
                        nodesinelem[index] = cur - vstart; 
                        index++;
                    }
                }
                
                DMPlexRestoreTransitiveClosure(mypetscmesh, i, PETSC_TRUE, &numpoints, &points);
                
                // Bring to out element node ordering:
                reordernodes(elemtypenum, nodesinelem);
                
                // Add element:
                int elementindexincurrenttype = myelements.add(elemtypenum, curvatureorder, nodesinelem);
                
                // Get the physical region and add the element:
                physicalregion* currentphysicalregion = myphysicalregions.get(curlabelvalue);
                currentphysicalregion->addelement(elemtypenum, elementindexincurrenttype);
            }
        }
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
    PetscViewerFormat format = PETSC_VIEWER_ASCII_VTK;
    PetscViewerASCIIOpen(PETSC_COMM_SELF, "petscview.txt", &v);
    PetscViewerPushFormat(v, format);
    DMView(mypetscmesh, v);
    PetscViewerDestroy(&v);
}

