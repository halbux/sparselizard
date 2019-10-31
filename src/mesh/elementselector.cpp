#include "elementselector.h"


void elementselector::prepare(bool isorientationdependent)
{
    // Sort the elements according to their total orientation.
    // Do it only if orientation dependent otherwise the disjoint
    // regions will not be sorted according to the order defined in
    // 'disjointregionnumbers'.
    if (isorientationdependent == true)
    {
        std::vector<int> renumberingvector;
        myalgorithm::stablesort(totalorientations, renumberingvector);

        std::vector<int> totalorientationsbackup = totalorientations;
        std::vector<int> disjointregionsbackup = disjointregions;
        std::vector<int> elemsbackup = elems;
        std::vector<int> originalindexesbackup = originalindexes;
        
        for (int i = 0; i < totalorientations.size(); i++)
        {
            totalorientations[i] = totalorientationsbackup[renumberingvector[i]];
            disjointregions[i] = disjointregionsbackup[renumberingvector[i]];
            elems[i] = elemsbackup[renumberingvector[i]];
            originalindexes[i] = originalindexesbackup[renumberingvector[i]];
        }
    }
    
    // Set the range begin and end:
    currentrangebegin = 0; currentrangeend = 0;
    currenttotalorientation = totalorientations[currentrangebegin];
    
    if (isorientationdependent == false)
        currentrangeend = totalorientations.size() - 1;
    else
    {
        while (currentrangeend < totalorientations.size() && totalorientations[currentrangeend] == currenttotalorientation)
            currentrangeend++;
        currentrangeend--;
    }   
}    

elementselector::elementselector(std::vector<int> disjointregionnumbers, bool isorientationdependent)
{
    // Have a coherent disjoint region ordering:
    std::sort(disjointregionnumbers.begin(), disjointregionnumbers.end());
    
    mydisjointregionnumbers = disjointregionnumbers;

    // Get the total number of elements in all disjoint regions for preallocation.
    int totalnumberofelements = 0;
    for (int i = 0; i < mydisjointregionnumbers.size(); i++)
        totalnumberofelements += universe::mymesh->getdisjointregions()->countelements(mydisjointregionnumbers[i]);
    
    // Create 'disjointregions', 'totalorientations' and 'elems':
    disjointregions.resize(totalnumberofelements);
    totalorientations.resize(totalnumberofelements);
    elems.resize(totalnumberofelements);
    originalindexes.resize(totalnumberofelements);
    std::iota(originalindexes.begin(), originalindexes.end(), 0);
    
    int currentindex = 0;
    elements* myelements = universe::mymesh->getelements();
    for (int i = 0; i < mydisjointregionnumbers.size(); i++)
    {
        int numelemindisjreg = universe::mymesh->getdisjointregions()->countelements(mydisjointregionnumbers[i]);
        int myelementtypenumber = universe::mymesh->getdisjointregions()->getelementtypenumber(mydisjointregionnumbers[i]);
        int rangebegin = universe::mymesh->getdisjointregions()->getrangebegin(mydisjointregionnumbers[i]);
        for (int j = 0; j < numelemindisjreg; j++)
        {
            disjointregions[currentindex+j] = mydisjointregionnumbers[i];
            totalorientations[currentindex+j] = myelements->gettotalorientation(myelementtypenumber, rangebegin+j);
            elems[currentindex+j] = rangebegin+j;
        }
        currentindex += numelemindisjreg;
    }
    
    prepare(isorientationdependent);
}    

elementselector::elementselector(std::vector<int> disjointregionnumbers, std::vector<int>& elemnums, bool isorientationdependent)
{
    // Have a coherent disjoint region ordering:
    std::sort(disjointregionnumbers.begin(), disjointregionnumbers.end());
    
    mydisjointregionnumbers = disjointregionnumbers;

    // Create 'disjointregions', 'totalorientations' and 'elems':
    disjointregions.resize(elemnums.size());
    totalorientations.resize(elemnums.size());
    elems.resize(elemnums.size());
    originalindexes.resize(elemnums.size());
    std::iota(originalindexes.begin(), originalindexes.end(), 0);
    
    int myelementtypenumber = universe::mymesh->getdisjointregions()->getelementtypenumber(mydisjointregionnumbers[0]);
    
    elements* myelements = universe::mymesh->getelements();
    for (int i = 0; i < elemnums.size(); i++)
    {
        disjointregions[i] = myelements->getdisjointregion(myelementtypenumber, elemnums[i]);
        totalorientations[i] = myelements->gettotalorientation(myelementtypenumber, elemnums[i]);
        elems[i] = elemnums[i];
    }
    
    prepare(isorientationdependent);
}

int elementselector::getelementdimension(void)
{
    return universe::mymesh->getdisjointregions()->getelementdimension(mydisjointregionnumbers[0]);
}

int elementselector::getelementtypenumber(void)
{
    return universe::mymesh->getdisjointregions()->getelementtypenumber(mydisjointregionnumbers[0]);
}

bool elementselector::next(void)
{
    currentrangebegin = currentrangeend+1;
    currentrangeend++;
    
    if (currentrangebegin >= totalorientations.size())
        return false;
    
    currenttotalorientation = totalorientations[currentrangebegin];
    
    while (currentrangeend < totalorientations.size() && totalorientations[currentrangeend] == currenttotalorientation)
        currentrangeend++;
    currentrangeend--;
    
    return true;
}    

void elementselector::selectallelements(void)
{
    // Set the range begin and end:
    currentrangebegin = 0; currentrangeend = totalorientations.size() - 1;
    currenttotalorientation = totalorientations[currentrangebegin];
    
    selecteddisjointregions = {};
}

void elementselector::selectdisjointregions(std::vector<int> disjregs) 
{ 
    std::sort(disjregs.begin(), disjregs.end()); 
    selecteddisjointregions = disjregs; 
}

int elementselector::countinselection(void)
{
    if (currentrangebegin >= totalorientations.size())
        return 0;
    else
    {
        // If all disjoint regions are requested:
        if (selecteddisjointregions.size() == 0)
            return currentrangeend - currentrangebegin + 1;
        else
        {
            // 'isdisjregrequested[disjreg]' is true if disjreg is in 'selecteddisjointregions'.
            std::vector<bool> isdisjregrequested(*std::max_element(mydisjointregionnumbers.begin(), mydisjointregionnumbers.end())+1, false);
            for (int i = 0; i < selecteddisjointregions.size(); i++)
                isdisjregrequested[selecteddisjointregions[i]] = true;
            
            int numinselection = 0;
            for (int i = currentrangebegin; i <= currentrangeend; i++)
            {
                if (isdisjregrequested[disjointregions[i]])
                    numinselection++;
            }
            return numinselection;
        }
    }
}    

std::vector<int> elementselector::getelementnumbers(void)
{
    int numinselection = countinselection();
    std::vector<int> elementnumbers(numinselection);
    
    // If all disjoint regions are requested:
    if (selecteddisjointregions.size() == 0)
    {
        for (int i = 0; i < countincurrentorientation(); i++)
            elementnumbers[i] = elems[currentrangebegin+i];
    }
    else
    {
        // 'isdisjregrequested[disjreg]' is true if disjreg is in 'disjregs'.
        std::vector<bool> isdisjregrequested(*std::max_element(mydisjointregionnumbers.begin(), mydisjointregionnumbers.end())+1, false);
        for (int i = 0; i < selecteddisjointregions.size(); i++)
            isdisjregrequested[selecteddisjointregions[i]] = true;
        
        int index = 0;
        for (int i = 0; i < countincurrentorientation(); i++)
        {
            if (isdisjregrequested[disjointregions[currentrangebegin+i]])
            {
                elementnumbers[index] = elems[currentrangebegin+i];
                index++;
            }
        }
    }
    return elementnumbers;
}    

std::vector<int> elementselector::getelementindexes(void)
{
    int numinselection = countinselection();
    std::vector<int> elementindexes(numinselection);
    
    // If all disjoint regions are requested:
    if (selecteddisjointregions.size() == 0)
    {
        for (int i = 0; i < countincurrentorientation(); i++)
            elementindexes[i] = i;
    }
    else
    {
        // 'isdisjregrequested[disjreg]' is true if disjreg is in 'disjregs'.
        std::vector<bool> isdisjregrequested(*std::max_element(mydisjointregionnumbers.begin(), mydisjointregionnumbers.end())+1, false);
        for (int i = 0; i < selecteddisjointregions.size(); i++)
            isdisjregrequested[selecteddisjointregions[i]] = true;
        
        int index = 0;
        for (int i = 0; i < countincurrentorientation(); i++)
        {
            if (isdisjregrequested[disjointregions[currentrangebegin+i]])
            {
                elementindexes[index] = i;
                index++;
            }
        }
    }
    return elementindexes;
}    

std::vector<int> elementselector::getoriginalindexes(void)
{
    int numinselection = countinselection();
    std::vector<int> origindexes(numinselection);
    
    // If all disjoint regions are requested:
    if (selecteddisjointregions.size() == 0)
    {
        for (int i = 0; i < countincurrentorientation(); i++)
            origindexes[i] = originalindexes[currentrangebegin+i];
    }
    else
    {
        // 'isdisjregrequested[disjreg]' is true if disjreg is in 'disjregs'.
        std::vector<bool> isdisjregrequested(*std::max_element(mydisjointregionnumbers.begin(), mydisjointregionnumbers.end())+1, false);
        for (int i = 0; i < selecteddisjointregions.size(); i++)
            isdisjregrequested[selecteddisjointregions[i]] = true;
        
        int index = 0;
        for (int i = 0; i < countincurrentorientation(); i++)
        {
            if (isdisjregrequested[disjointregions[currentrangebegin+i]])
            {
                origindexes[index] = originalindexes[currentrangebegin+i];
                index++;
            }
        }
    }
    return origindexes;
}    

elementselector elementselector::extractselection(void)
{
    elementselector extracted;
    
    if (selecteddisjointregions.size() == 0)
        extracted.mydisjointregionnumbers = mydisjointregionnumbers;
    else
        extracted.mydisjointregionnumbers = selecteddisjointregions;
    extracted.currenttotalorientation = currenttotalorientation;
    
    // There is only a single orientation in 'extracted':
    int numinselection = countinselection();
    extracted.currentrangebegin = 0;
    extracted.currentrangeend = numinselection-1;
    
    // Create the three containers in 'extracted':
    extracted.elems.resize(numinselection);
    extracted.totalorientations.resize(numinselection);
    extracted.disjointregions.resize(numinselection);
    extracted.originalindexes.resize(numinselection);
    
    std::vector<int> selectedindexes = getelementindexes();
    for (int i = 0; i < numinselection; i++)
    {
        extracted.elems[i] = elems[currentrangebegin+selectedindexes[i]];
        extracted.totalorientations[i] = totalorientations[currentrangebegin+selectedindexes[i]];
        extracted.disjointregions[i] = disjointregions[currentrangebegin+selectedindexes[i]];
        extracted.originalindexes[i] = originalindexes[currentrangebegin+selectedindexes[i]];
    }
    
    return extracted;
}


