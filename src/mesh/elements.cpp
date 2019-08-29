#include "elements.h"
#include "geotools.h"


elements::elements(nodes& inputnodes, physicalregions& inputphysicalregions, disjointregions& inputdisjointregions)
{
    mynodes = &inputnodes;
    myphysicalregions = &inputphysicalregions;
    mydisjointregions = &inputdisjointregions;
}

int elements::add(int elementtypenumber, int curvatureorder, std::vector<int>& nodelist)
{
    // Point elements have a number equal to the node defining them
    // and should not be added to 'subelementsinelements'.
    if (elementtypenumber == 0)
        return nodelist[0];

    // Define the curvature order for the first time:
    if (mycurvatureorder == -1)
    {
        mycurvatureorder = curvatureorder;
        // Define as well the following quantities now that the curvature order is known:
        for (int i = 0; i <= 7; i++)
        {
            element myelement(i,curvatureorder);
            
            numberofsubelementsineveryelement[i][0] = myelement.countcurvednodes();
            numberofsubelementsineveryelement[i][1] = myelement.countedges();
            numberofsubelementsineveryelement[i][2] = myelement.counttriangularfaces();
            numberofsubelementsineveryelement[i][3] = myelement.countquadrangularfaces();
        }
    }
    if (mycurvatureorder != curvatureorder)
    {
        std::cout << "Error in 'elements' object: the mesh can contain only a single curvature order." << std::endl;
        abort();
    }
        
    subelementsinelements[elementtypenumber][0].insert(subelementsinelements[elementtypenumber][0].end(), nodelist.begin(), nodelist.end());
    return subelementsinelements[elementtypenumber][0].size()/nodelist.size() - 1;
}

void elements::shift(double xshift, double yshift, double zshift)
{
    for (int i = 0; i < barycenters.size(); i++)
    {
        for (int j = 0; j < barycenters[i].size()/3; j++)
        {
            barycenters[i][3*j+0] += xshift;
            barycenters[i][3*j+1] += yshift;
            barycenters[i][3*j+2] += zshift;
        }
    }
}

void elements::rotate(double alphax, double alphay, double alphaz)
{
    for (int i = 0; i < barycenters.size(); i++)
    {
        if (barycenters[i].size() > 0)
            geotools::rotate(alphax, alphay, alphaz, &barycenters[i]);
    }
}

void elements::scale(double xscale, double yscale, double zscale)
{
    for (int i = 0; i < barycenters.size(); i++)
    {
        for (int j = 0; j < barycenters[i].size()/3; j++)
        {
            barycenters[i][3*j+0] *= xscale;
            barycenters[i][3*j+1] *= yscale;
            barycenters[i][3*j+2] *= zscale;
        }
    }
}

int elements::getsubelement(int subelementtypenumber, int elementtypenumber, int elementnumber, int subelementindex)
{
    if (elementtypenumber == subelementtypenumber)
        return elementnumber;
    else
        return subelementsinelements[elementtypenumber][subelementtypenumber][elementnumber*numberofsubelementsineveryelement[elementtypenumber][subelementtypenumber]+subelementindex];
}

int elements::getdisjointregion(int elementtypenumber, int elementnumber)
{
    return indisjointregion[elementtypenumber][elementnumber];
}

int elements::gettotalorientation(int elementtypenumber, int elementnumber)
{
    if (elementtypenumber != 0)
        return totalorientations[elementtypenumber][elementnumber];
    else
        return 0;
}


int elements::count(int elementtypenumber)
{
    if (elementtypenumber == 0)
        return mynodes->count();
    
    if (subelementsinelements[elementtypenumber][0].size() == 0)
        return 0;
    
    return subelementsinelements[elementtypenumber][0].size()/numberofsubelementsineveryelement[elementtypenumber][0];
}

int elements::getcurvatureorder(void)
{
    return mycurvatureorder;
}

void elements::populateedgesatnodes(void)
{
    // Get the number of nodes:
    int numnodes = count(0);
    // Get the number of edges:
    int numedges = count(1);
    // Get the number of curved nodes per edge:
    int numcurvednodes = subelementsinelements[1][0].size()/numedges;


    // Preallocate vectors:
    adressedgesatnodes = std::vector<int>(numnodes,0);
    edgesatnodes = std::vector<int>(numedges*2);


    // Vector to count the number of edges touching every node:
    std::vector<int> numedgesonnode(numnodes,0);

    // Loop on all edges:
    for (int e = 0; e < numedges; e++)
    {
        // Corner nodes in current edge:
        int edgefirstnode = subelementsinelements[1][0][e*numcurvednodes+0];
        int edgelastnode = subelementsinelements[1][0][e*numcurvednodes+1];

        numedgesonnode[edgefirstnode]++;
        numedgesonnode[edgelastnode]++;
    }

    for (int n = 1; n < numnodes; n++)
        adressedgesatnodes[n] = adressedgesatnodes[n-1] + numedgesonnode[n-1];

    std::vector<int> currentedgeinnode(numnodes,0);    
    for (int e = 0; e < numedges; e++)
    {
        // Corner nodes in current edge:
        int edgefirstnode = subelementsinelements[1][0][e*numcurvednodes+0];
        int edgelastnode = subelementsinelements[1][0][e*numcurvednodes+1];

        edgesatnodes[adressedgesatnodes[edgefirstnode]+currentedgeinnode[edgefirstnode]] = e;
        currentedgeinnode[edgefirstnode]++;
        edgesatnodes[adressedgesatnodes[edgelastnode]+currentedgeinnode[edgelastnode]] = e;
        currentedgeinnode[edgelastnode]++;
    }
}

int elements::countedgesonnode(int nodenumber)
{
    if (edgesatnodes.size() == 0)
        populateedgesatnodes();

    if (nodenumber+1 < adressedgesatnodes.size())
        return adressedgesatnodes[nodenumber+1]-adressedgesatnodes[nodenumber];
    else
        return edgesatnodes.size()-adressedgesatnodes[nodenumber];
}

std::vector<int> elements::getedgesonnode(int nodenumber)
{
    int numedgesatnode = countedgesonnode(nodenumber);

    std::vector<int> output(numedgesatnode);
    for (int i = 0; i < numedgesatnode; i++)
        output[i] = edgesatnodes[adressedgesatnodes[nodenumber]+i];

    return output;
}

std::vector<double> elements::getnodecoordinates(int elementtypenumber, int elementnumber, int xyz)
{
    std::vector<double>* nodecoordinates = mynodes->getcoordinates();
    
    if (elementtypenumber == 0)
        return {nodecoordinates->at(3*elementnumber+xyz)};
    
    element myelement(elementtypenumber, mycurvatureorder);
    int curvednumberofnodes = myelement.countcurvednodes();
    
    std::vector<double> nodecoords(curvednumberofnodes);

    for (int node = 0; node < curvednumberofnodes; node++)
        nodecoords[node] = nodecoordinates->at(3*subelementsinelements[elementtypenumber][0][elementnumber*curvednumberofnodes+node]+xyz);

    return nodecoords;
}

std::vector<double>* elements::getbarycenters(int elementtypenumber)
{
    // If not yet populated for the element type:
    if (barycenters[elementtypenumber].size() == 0)
        barycenters[elementtypenumber] = computebarycenters(elementtypenumber);
    
    return &(barycenters[elementtypenumber]);
}

std::vector<double>* elements::getsphereradius(int elementtypenumber)
{
    std::vector<double>* mybarys = getbarycenters(elementtypenumber);

    // If not yet populated for the element type:
    if (sphereradius[elementtypenumber].size() == 0)
    {
        sphereradius[elementtypenumber].resize(count(elementtypenumber));
    
        for (int i = 0; i < count(elementtypenumber); i++)
        {
            double maxdist = 0;
        
            std::vector<double> xnodes = getnodecoordinates(elementtypenumber, i, 0);
            std::vector<double> ynodes = getnodecoordinates(elementtypenumber, i, 1);
            std::vector<double> znodes = getnodecoordinates(elementtypenumber, i, 2);
            
            for (int j = 0; j < xnodes.size(); j++)
            {
                double curdist = std::sqrt( std::pow(mybarys->at(3*i+0)-xnodes[j], 2) + std::pow(mybarys->at(3*i+1)-ynodes[j], 2) + std::pow(mybarys->at(3*i+2)-znodes[j], 2) );
                if (curdist > maxdist)
                    maxdist = curdist;
            }
            
            sphereradius[elementtypenumber][i] = maxdist;
        }
    }
    
    return &(sphereradius[elementtypenumber]);
}

void elements::printnumber(void)
{
    std::cout << std::endl;
    for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
    {
        element myelement(elementtypenumber);
        std::cout << "Number of " << std::left << std::setw(12) << myelement.gettypenameconjugation(2) << " " << count(elementtypenumber) << std::endl;
    }
    std::cout << std::endl;
}

void elements::printsubelements(void)
{
    for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
    {
        element myelement(elementtypenumber, mycurvatureorder);

        if (subelementsinelements[elementtypenumber][0].size() != 0)
        {
            std::cout << "Points in " << myelement.gettypenameconjugation(2) << ":" << std::endl;
            for (int i = 0; i < subelementsinelements[elementtypenumber][0].size()/myelement.countcurvednodes(); i++)
            {
                std::cout << std::left << std::setw(8) << i << " ";
                for (int j = 0; j < myelement.countcurvednodes(); j++)
                    std::cout << std::left << std::setw(6) << subelementsinelements[elementtypenumber][0][i*myelement.countcurvednodes()+j] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        if (subelementsinelements[elementtypenumber][1].size() != 0)
        {
            std::cout << "Lines in " << myelement.gettypenameconjugation(2) << ":" << std::endl;
            for (int i = 0; i < subelementsinelements[elementtypenumber][1].size()/myelement.countedges(); i++)
            {
                std::cout << std::left << std::setw(8) << i << " ";
                for (int j = 0; j < myelement.countedges(); j++)
                    std::cout << std::left << std::setw(6) << subelementsinelements[elementtypenumber][1][i*myelement.countedges()+j] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        if (subelementsinelements[elementtypenumber][2].size() != 0)
        {
            std::cout << "Triangles in " << myelement.gettypenameconjugation(2) << ":" << std::endl;
            for (int i = 0; i < subelementsinelements[elementtypenumber][2].size()/myelement.counttriangularfaces(); i++)
            {
                std::cout << std::left << std::setw(8) << i << " ";
                for (int j = 0; j < myelement.counttriangularfaces(); j++)
                    std::cout << std::left << std::setw(6) << subelementsinelements[elementtypenumber][2][i*myelement.counttriangularfaces()+j] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        if (subelementsinelements[elementtypenumber][3].size() != 0)
        {
            std::cout << "Quadrangles in " << myelement.gettypenameconjugation(2) << ":" << std::endl;
            for (int i = 0; i < subelementsinelements[elementtypenumber][3].size()/myelement.countquadrangularfaces(); i++)
            {
                std::cout << std::left << std::setw(8) << i << " ";
                for (int j = 0; j < myelement.countquadrangularfaces(); j++)
                    std::cout << std::left << std::setw(6) << subelementsinelements[elementtypenumber][3][i*myelement.countquadrangularfaces()+j] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
}

void elements::printtotalorientations(void)
{
    std::cout << std::endl;
    for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
    {
        if (totalorientations[elementtypenumber].size() == 0)
            continue;
        
        element myelement(elementtypenumber);
        std::cout << "Total orientations for all " << myelement.gettypenameconjugation(2) << ":" << std::endl;
        
        for (int i = 0; i < totalorientations[elementtypenumber].size(); i++)
            std::cout << std::left << std::setw(12) << i << std::left << std::setw(12) << totalorientations[elementtypenumber][i] << std::endl;
    }
    std::cout << std::endl;
}

std::vector<double> elements::computebarycenters(int elementtypenumber)
{
    std::vector<double>* nodecoordinates = mynodes->getcoordinates();

    // The barycenter of point elements is the nodes coordinates:
    if (elementtypenumber == 0)
        return *nodecoordinates;
        
    element myelement(elementtypenumber, mycurvatureorder);
    int curvednumberofnodes = myelement.countcurvednodes();
    
    // Preallocate the barycenter coordinates vector:
    std::vector<double> barycentercoordinates(3 * count(elementtypenumber),0);
    // Compute the barycenters on the straight element:
    for (int elem = 0; elem < count(elementtypenumber); elem++)
    {
        for (int node = 0; node < myelement.countnodes(); node++)
        {
            barycentercoordinates[3*elem+0] += nodecoordinates->at(3*subelementsinelements[elementtypenumber][0][elem*curvednumberofnodes+node]+0);
            barycentercoordinates[3*elem+1] += nodecoordinates->at(3*subelementsinelements[elementtypenumber][0][elem*curvednumberofnodes+node]+1);
            barycentercoordinates[3*elem+2] += nodecoordinates->at(3*subelementsinelements[elementtypenumber][0][elem*curvednumberofnodes+node]+2);
        }
        barycentercoordinates[3*elem+0] = barycentercoordinates[3*elem+0]/myelement.countnodes();
        barycentercoordinates[3*elem+1] = barycentercoordinates[3*elem+1]/myelement.countnodes();
        barycentercoordinates[3*elem+2] = barycentercoordinates[3*elem+2]/myelement.countnodes();
    }
    
    return barycentercoordinates;
}

std::vector<int> elements::sortbybarycenters(int elementtypenumber)
{
    // Get the barycenters of the straight version of all elements of type 'elementtypenumber':
    std::vector<double> barycentercoordinates = computebarycenters(elementtypenumber);
    
    // 'reorderingvector' gives the relation between the indexes before and after element sorting:
    std::vector<int> reorderingvector;
    myalgorithm::stablecoordinatesort(mynodes->getnoisethreshold(), barycentercoordinates, reorderingvector);

    // sortedcoordinates = coordinates(reorderingvector,:).
    // sortedcoordinates(renumberingvector,:) = coordinates.
    std::vector<int> renumberingvector(count(elementtypenumber));
    for (int i = 0; i < count(elementtypenumber); i++)
        renumberingvector[reorderingvector[i]] = i;
    
    // Renumber and reorder the elements in all containers:
    renumber(elementtypenumber, renumberingvector);
    reorder(elementtypenumber, reorderingvector);
    
    return renumberingvector;
}

std::vector<int> elements::removeduplicates(int elementtypenumber)
{
    // For point elements (i.e. to remove node duplicates):
    if (elementtypenumber == 0)
    {
        std::vector<int> renumberingvector = mynodes->removeduplicates();
        renumber(0, renumberingvector);
        return renumberingvector;
    }
    
    int numberofcurvednodes = numberofsubelementsineveryelement[elementtypenumber][0];
    int numberoflines = numberofsubelementsineveryelement[elementtypenumber][1];
    int numberoftriangles = numberofsubelementsineveryelement[elementtypenumber][2];
    int numberofquadrangles = numberofsubelementsineveryelement[elementtypenumber][3];
    
    std::vector<double> noisethreshold = mynodes->getnoisethreshold();
    std::vector<double> barycentercoordinates = computebarycenters(elementtypenumber);
    
    // 'elementrenumbering' will give the renumbering corresponding to removed duplicates:
    std::vector<int> elementrenumbering;
    int numberofnonduplicates = myalgorithm::removeduplicatedcoordinates(noisethreshold, barycentercoordinates, elementrenumbering);
    
    for (int i = 0; i < elementrenumbering.size(); i++)
    {
        if (elementrenumbering[i] != i)
        {
            for (int j = 0; j < numberofcurvednodes; j++)
                subelementsinelements[elementtypenumber][0][elementrenumbering[i]*numberofcurvednodes+j] = subelementsinelements[elementtypenumber][0][i*numberofcurvednodes+j];
            // Watch out: 'subelementsinelements[1][1]' does not contain anything!
            if (subelementsinelements[elementtypenumber][1].size() != 0)
            {
                for (int j = 0; j < numberoflines; j++)
                    subelementsinelements[elementtypenumber][1][elementrenumbering[i]*numberoflines+j] = subelementsinelements[elementtypenumber][1][i*numberoflines+j];
            }
            if (subelementsinelements[elementtypenumber][2].size() != 0)
            {
                for (int j = 0; j < numberoftriangles; j++)
                    subelementsinelements[elementtypenumber][2][elementrenumbering[i]*numberoftriangles+j] = subelementsinelements[elementtypenumber][2][i*numberoftriangles+j];
            }
            if (subelementsinelements[elementtypenumber][3].size() != 0)
            {
                for (int j = 0; j < numberofquadrangles; j++)
                    subelementsinelements[elementtypenumber][3][elementrenumbering[i]*numberofquadrangles+j] = subelementsinelements[elementtypenumber][3][i*numberofquadrangles+j];
            }
        }
    }
    // Shrink to fit:
    subelementsinelements[elementtypenumber][0].resize(numberofnonduplicates*numberofcurvednodes);
    if (subelementsinelements[elementtypenumber][1].size() != 0)
        subelementsinelements[elementtypenumber][1].resize(numberofnonduplicates*numberoflines);
    if (subelementsinelements[elementtypenumber][2].size() != 0)
        subelementsinelements[elementtypenumber][2].resize(numberofnonduplicates*numberoftriangles);
    if (subelementsinelements[elementtypenumber][3].size() != 0)
        subelementsinelements[elementtypenumber][3].resize(numberofnonduplicates*numberofquadrangles);

    renumber(elementtypenumber, elementrenumbering);
    
    return elementrenumbering;
}

void elements::renumber(int elementtypenumber, std::vector<int>& renumberingvector)
{   
    for (int typenum = 0; typenum <= 7; typenum++)
    {
        switch (elementtypenumber)
        {
            // Renumber the point elements (i.e. the nodes):
            case 0:
                for (int i = 0; i < subelementsinelements[typenum][0].size(); i++)
                    subelementsinelements[typenum][0][i] = renumberingvector[subelementsinelements[typenum][0][i]];
                break;
            // Renumber the line elements:
            case 1:
                for (int i = 0; i < subelementsinelements[typenum][1].size(); i++)
                    subelementsinelements[typenum][1][i] = renumberingvector[subelementsinelements[typenum][1][i]];
                break;
            // Renumber the triangle elements:
            case 2:
                for (int i = 0; i < subelementsinelements[typenum][2].size(); i++)
                    subelementsinelements[typenum][2][i] = renumberingvector[subelementsinelements[typenum][2][i]];
                break;
            // Renumber the quadrangle elements:
            case 3:
                for (int i = 0; i < subelementsinelements[typenum][3].size(); i++)
                    subelementsinelements[typenum][3][i] = renumberingvector[subelementsinelements[typenum][3][i]];
                break;
        }
    }
}

void elements::reorder(int elementtypenumber, std::vector<int> &elementreordering)
{
    // For point elements (i.e. nodes):
    if (elementtypenumber == 0)
    {
        mynodes->reorder(elementreordering);
        return;
    }
    
    int numberofnodes = numberofsubelementsineveryelement[elementtypenumber][0];
    int numberoflines = numberofsubelementsineveryelement[elementtypenumber][1];
    int numberoftriangles = numberofsubelementsineveryelement[elementtypenumber][2];
    int numberofquadrangles = numberofsubelementsineveryelement[elementtypenumber][3];
    
    std::vector<int> pointsinelementspart = subelementsinelements[elementtypenumber][0];
    for (int i = 0; i < subelementsinelements[elementtypenumber][0].size()/numberofnodes; i++)
    {
        for (int j = 0; j < numberofnodes; j++)
            subelementsinelements[elementtypenumber][0][numberofnodes*i+j] = pointsinelementspart[numberofnodes*elementreordering[i]+j]; 
    }
    std::vector<int> linesinelementspart = subelementsinelements[elementtypenumber][1];
    for (int i = 0; i < subelementsinelements[elementtypenumber][1].size()/numberoflines; i++)
    {
        for (int j = 0; j < numberoflines; j++)
            subelementsinelements[elementtypenumber][1][numberoflines*i+j] = linesinelementspart[numberoflines*elementreordering[i]+j]; 
    }
    if (numberoftriangles != 0)
    {
        std::vector<int> trianglesinelementspart = subelementsinelements[elementtypenumber][2];
        for (int i = 0; i < subelementsinelements[elementtypenumber][2].size()/numberoftriangles; i++)
        {
            for (int j = 0; j < numberoftriangles; j++)
                subelementsinelements[elementtypenumber][2][numberoftriangles*i+j] = trianglesinelementspart[numberoftriangles*elementreordering[i]+j]; 
        }
    }
    if (numberofquadrangles != 0)
    {
        std::vector<int> quadranglesinelementspart = subelementsinelements[elementtypenumber][3];
        for (int i = 0; i < subelementsinelements[elementtypenumber][3].size()/numberofquadrangles; i++)
        {
            for (int j = 0; j < numberofquadrangles; j++)
                subelementsinelements[elementtypenumber][3][numberofquadrangles*i+j] = quadranglesinelementspart[numberofquadrangles*elementreordering[i]+j]; 
        }
    }
}

void elements::explode(void)
{                            
    // Add all new elements. Loop on all elements with increasing dimension 
    // to avoid defining too many duplicates.
    // Skip lines (type 1) since there is nothing to add for lines.
    for (int elementtypenumber = 2; elementtypenumber <= 7; elementtypenumber++)
    {
        element myelement(elementtypenumber, mycurvatureorder);
        int curvednumberofnodes = myelement.countcurvednodes();
        
        std::vector<int> facesdefinitionsbasedonedges = myelement.getfacesdefinitionsbasedonedges();

        // Skip if there are no elements of the current type:
        if (count(elementtypenumber) == 0)
            continue;

        std::vector<int> currentnodes(curvednumberofnodes,0);

        // Loop on all elements of the current type:
        for (int elem = 0; elem < count(elementtypenumber); elem++)
        {
            // Define the nodes of the element object:
            for (int i = 0; i < curvednumberofnodes; i++)
                currentnodes[i] = subelementsinelements[elementtypenumber][0][elem*curvednumberofnodes+i];
            myelement.setnodes(currentnodes);

            // Add all line subelements - only for 2D and 3D elements:
            int firstnewlinenumber = count(1);
            for (int line = 0; line < myelement.countedges(); line++)
            {
                std::vector<int> nodesinline = myelement.getnodesinline(line);
                int newlinenumber = add(1, mycurvatureorder, nodesinline);
                subelementsinelements[elementtypenumber][1].push_back(newlinenumber);
            }
            // Add all subsurfaces - only for 3D elements:
            if (myelement.getelementdimension() == 3)
            {
                for (int triangle = 0; triangle < myelement.counttriangularfaces(); triangle++)
                {
                    std::vector<int> nodesintriangle = myelement.getnodesintriangle(triangle);
                    int newtrianglenumber = add(2, mycurvatureorder, nodesintriangle);
                    subelementsinelements[elementtypenumber][2].push_back(newtrianglenumber);
                    // Also link the triangle to its lines. All lines have already been defined before above.
                    for (int triangleline = 0; triangleline < 3; triangleline++)
                        subelementsinelements[2][1].push_back(firstnewlinenumber + std::abs(facesdefinitionsbasedonedges[3*triangle+triangleline])-1);
                }    
                for (int quadrangle = 0; quadrangle < myelement.countquadrangularfaces(); quadrangle++)
                {
                    std::vector<int> nodesinquadrangle = myelement.getnodesinquadrangle(quadrangle);
                    int newquadranglenumber = add(3, mycurvatureorder, nodesinquadrangle);
                    subelementsinelements[elementtypenumber][3].push_back(newquadranglenumber);
                    // Also link the quadrangle to its lines (defined before):
                    for (int quadrangleline = 0; quadrangleline < 4; quadrangleline++)
                        subelementsinelements[3][1].push_back(firstnewlinenumber + std::abs(facesdefinitionsbasedonedges[3*myelement.counttriangularfaces()+4*quadrangle+quadrangleline])-1);
                }    
            }
        }
    }
}

void elements::definedisjointregions(void)
{
    
    int numberofphysicalregions = myphysicalregions->count();
    
    // Define 'isinphysicalregion[elementtypenumber]' to hold a vector
    // whose (elem*numphysreg+i)th entry is true if the element 'elem' 
    // is included in the ith physical region.
    std::vector<std::vector<bool>> isinphysicalregion(8);
    for (int typenum = 0; typenum <= 7; typenum++)
        isinphysicalregion[typenum] = std::vector<bool>(count(typenum)*numberofphysicalregions);
    
    // Write the physical regions to the elements. 
    for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
    {
        physicalregion* currentphysicalregion = myphysicalregions->getatindex(physregindex);
        std::vector<std::vector<int>>* elementsinphysicalregion = currentphysicalregion->getelementlist();
        
        for (int typenum = 0; typenum <= 7; typenum++)
        {
            // Iterate on all elements of the given type:
            for (int i = 0; i < (*elementsinphysicalregion)[typenum].size(); i++)
                isinphysicalregion[typenum][(*elementsinphysicalregion)[typenum][i] * numberofphysicalregions + physregindex] = true;
        }
    }
        
    
    // Propagate the physical region memberships from elements to the subelements.
    for (int typenum = 1; typenum <= 7; typenum++)
    {
        int numberofcurvednodes = numberofsubelementsineveryelement[typenum][0];
        int numberoflines = numberofsubelementsineveryelement[typenum][1];
        int numberoftriangles = numberofsubelementsineveryelement[typenum][2];
        int numberofquadrangles = numberofsubelementsineveryelement[typenum][3];

        for (int elem = 0; elem < subelementsinelements[typenum][0].size()/numberofcurvednodes; elem++)
        {
            for (int i = 0; i < numberofcurvednodes; i++)
            {
                int subelem = subelementsinelements[typenum][0][elem*numberofcurvednodes+i];
                for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                {
                    if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                        isinphysicalregion[0][subelem*numberofphysicalregions+physregindex] = true;
                }
            }
        }
        
        for (int elem = 0; elem < subelementsinelements[typenum][1].size()/numberoflines; elem++)
        {
            for (int i = 0; i < numberoflines; i++)
            {
                int subelem = subelementsinelements[typenum][1][elem*numberoflines+i];
                for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                {
                    if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                        isinphysicalregion[1][subelem*numberofphysicalregions+physregindex] = true;
                }
            }
        }

        if (numberoftriangles != 0)
        {
            for (int elem = 0; elem < subelementsinelements[typenum][2].size()/numberoftriangles; elem++)
            {
                for (int i = 0; i < numberoftriangles; i++)
                {
                    int subelem = subelementsinelements[typenum][2][elem*numberoftriangles+i];
                    for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                    {
                        if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                            isinphysicalregion[2][subelem*numberofphysicalregions+physregindex] = true;
                    }
                }
            }
        }
        
        if (numberofquadrangles != 0)
        {
            for (int elem = 0; elem < subelementsinelements[typenum][3].size()/numberofquadrangles; elem++)
            {
                for (int i = 0; i < numberofquadrangles; i++)
                {
                    int subelem = subelementsinelements[typenum][3][elem*numberofquadrangles+i];
                    for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                    {
                        if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                            isinphysicalregion[3][subelem*numberofphysicalregions+physregindex] = true;
                    }
                }
            }
        }
    }
        

    // Now define the disjoint regions based on 'isinphysicalregion':
    indisjointregion.resize(8);
    for (int typenum = 0; typenum <= 7; typenum++)
        indisjointregion[typenum].resize(count(typenum));
    // We need to know if a node is a corner node or not:
    std::vector<bool> isnodeacornernode = iscornernode();
    // Define the disjoint regions and fill in 'indisjointregion':
    for (int typenum = 0; typenum <= 7; typenum++)
    {
        if (count(typenum) == 0)
            continue;

        std::vector<bool> temp(numberofphysicalregions);    

        for (int elem = 0; elem < count(typenum); elem++)
        {   
            // Only corner nodes matter for the disjoint regions 
            // (no dof will ever be associated to curvature nodes).
            if (typenum != 0 || isnodeacornernode[elem])
            {
                for (int i = 0; i < numberofphysicalregions; i++)
                    temp[i] = isinphysicalregion[typenum][elem*numberofphysicalregions+i];

                indisjointregion[typenum][elem] = mydisjointregions->add(typenum, temp);
            }
            else
                indisjointregion[typenum][elem] = -1;
        }
    }
}

void elements::reorderbydisjointregions(void)
{
    for (int typenum = 0; typenum <= 7; typenum++)
    {
        std::vector<int> elementreordering;
        myalgorithm::stablesort(indisjointregion[typenum], elementreordering);
        
        std::vector<int> elementrenumbering(count(typenum));
        for (int i = 0; i < count(typenum); i++)
            elementrenumbering[elementreordering[i]] = i;
        
        reorder(typenum, elementreordering);
        renumber(typenum, elementrenumbering);

        // Reorder 'indisjointregion' as well:
        std::vector<int> indisjointregionpart = indisjointregion[typenum];
        for (int i = 0; i < indisjointregion[typenum].size(); i++)
            indisjointregion[typenum][i] = indisjointregionpart[elementreordering[i]]; 
    }
}

void elements::definedisjointregionsranges(void)
{
    for (int typenum = 0; typenum <= 7; typenum++)
    {
        for (int i = 0; i < indisjointregion[typenum].size(); i++)
        {
            if (i == 0)
                mydisjointregions->setrangebegin(indisjointregion[typenum][i], 0);
            if (i != 0 && indisjointregion[typenum][i] != indisjointregion[typenum][i-1])
            {
                mydisjointregions->setrangeend(indisjointregion[typenum][i-1], i-1);
                mydisjointregions->setrangebegin(indisjointregion[typenum][i], i);
            }
            if (i == indisjointregion[typenum].size()-1)
                mydisjointregions->setrangeend(indisjointregion[typenum][i], indisjointregion[typenum].size()-1);
        }
    }
}

std::vector<bool> elements::iscornernode(void)
{
    std::vector<bool> output(count(0), true);
    for (int typenum = 1; typenum <= 7; typenum++)
    {
        element curvedelement(typenum, getcurvatureorder());
        int numcurvednodes = curvedelement.countcurvednodes();
        int numcornernodes = curvedelement.countnodes();

        for (int i = 0; i < count(typenum); i++)
        {    
            // The first 'numcornernodes' are corner nodes, the remaining ones not.
            for (int node = numcornernodes; node < numcurvednodes; node++)
                output[subelementsinelements[typenum][0][i*numcurvednodes+node]] = false;
        }
    }
    return output;
}

void elements::tostandardorientation(void)
{
    // Only straight elements are supported for now:
    if (getcurvatureorder() != 1)
    {
        std::cout << "Note: curved elements are not reoriented for now (some functions may thus perform slower)" << std::endl;
        return;
    }

    // Loop on all element types except the point element (type 0):
    for (int elementtypenumber = 1; elementtypenumber <= 7; elementtypenumber++)
    {    
        if (subelementsinelements[elementtypenumber][0].size() == 0)
            continue;
        
        element myelement(elementtypenumber, mycurvatureorder);
        
        int numelemofcurrenttype = count(elementtypenumber);
        
        int numberofcurvednodes = myelement.countcurvednodes();
        int numberofedges = myelement.countedges();
        int numberoftriangularfaces = myelement.counttriangularfaces();
        int numberofquadrangularfaces = myelement.countquadrangularfaces();
        
        // Loop on all elements:
        std::vector<int> currentnodes(numberofcurvednodes);
        std::vector<int> currentedges(numberofedges);
        std::vector<int> currenttriangularfaces(numberoftriangularfaces);
        std::vector<int> currentquadrangularfaces(numberofquadrangularfaces);
        
        for (int i = 0; i < numelemofcurrenttype; i++)
        {
            for (int j = 0; j < numberofcurvednodes; j++)
                currentnodes[j] = subelementsinelements[elementtypenumber][0][i*numberofcurvednodes+j];
            myelement.setnodes(currentnodes);
            // This gives the corner nodes reordering:
            std::vector<int> nodereordering = myelement.getstandardorientationreordering();
            for (int j = 0; j < numberofcurvednodes; j++)
                subelementsinelements[elementtypenumber][0][i*numberofcurvednodes+j] = currentnodes[nodereordering[j]];
            // Reorder the edges:
            if (numberofedges > 0 && elementtypenumber != 1)
            {
                for (int j = 0; j < numberofedges; j++)
                    currentedges[j] = subelementsinelements[elementtypenumber][1][i*numberofedges+j];
                std::vector<int> edgereordering = myelement.getedgesreordering(nodereordering);
                for (int j = 0; j < numberofedges; j++)
                    subelementsinelements[elementtypenumber][1][i*numberofedges+j] = currentedges[edgereordering[j]];
            }
            // Reorder the triangular faces:
            if (numberoftriangularfaces > 0 && elementtypenumber != 2)
            {
                for (int j = 0; j < numberoftriangularfaces; j++)
                    currenttriangularfaces[j] = subelementsinelements[elementtypenumber][2][i*numberoftriangularfaces+j];
                std::vector<int> triangularfacereordering = myelement.gettriangularfacesreordering(nodereordering);
                for (int j = 0; j < numberoftriangularfaces; j++)
                    subelementsinelements[elementtypenumber][2][i*numberoftriangularfaces+j] = currenttriangularfaces[triangularfacereordering[j]];
            }
            // Reorder the quadrangular faces:
            if (numberofquadrangularfaces > 0 && elementtypenumber != 3)
            {
                for (int j = 0; j < numberofquadrangularfaces; j++)
                    currentquadrangularfaces[j] = subelementsinelements[elementtypenumber][3][i*numberofquadrangularfaces+j];
                std::vector<int> quadrangularfacereordering = myelement.getquadrangularfacesreordering(nodereordering);
                for (int j = 0; j < numberofquadrangularfaces; j++)
                    subelementsinelements[elementtypenumber][3][i*numberofquadrangularfaces+j] = currentquadrangularfaces[quadrangularfacereordering[j]];
            }
        }
    }
}

void elements::orient(void)
{
    // Loop on all element types except the point element (type 0):
    for (int elementtypenumber = 1; elementtypenumber <= 7; elementtypenumber++)
    {    
        if (subelementsinelements[elementtypenumber][0].size() == 0)
            continue;
        
        element myelement(elementtypenumber, mycurvatureorder);
        
        int numberofnodes = myelement.countnodes();
        int numberofcurvednodes = myelement.countcurvednodes();
        int numelemofcurrenttype = count(elementtypenumber);
        
        totalorientations[elementtypenumber].resize(numelemofcurrenttype);
    
        // Loop on all elements:
        std::vector<int> currentnodes(numberofnodes);
        for (int i = 0; i < numelemofcurrenttype; i++)
        {
            // Get only the corner nodes, not the curvature-related nodes:
            for (int j = 0; j < numberofnodes; j++)
                currentnodes[j] = subelementsinelements[elementtypenumber][0][i*numberofcurvednodes+j];
            
            totalorientations[elementtypenumber][i] = orientation::gettotalorientation(elementtypenumber, currentnodes); 
        }
    }
}
