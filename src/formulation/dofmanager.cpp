#include "dofmanager.h"


void dofmanager::addtostructure(shared_ptr<rawfield> fieldtoadd, std::vector<int> selecteddisjointregions)
{  
    disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();

    // Find the field index of 'fieldtoadd' (if present):
    int fieldindex = -1;
    for (int i = 0; i < myfields.size(); i++)
    {
        if (myfields[i].get() == fieldtoadd.get())
        {
            fieldindex = i;
            break;
        }
    }
    // Add the field to the structure if not existing.
    // Add an entry for every disjoint region at the same time.
    if (fieldindex == -1)
    {
        fieldindex = myfields.size();
        myfields.push_back(fieldtoadd);
        int numberofdisjointregions = mydisjointregions->count();
        std::vector<std::vector< int >> temp(numberofdisjointregions, std::vector< int >(0));
        rangebegin.push_back(temp);
        rangeend.push_back(temp);
    }
    
    // Add an entry for every form function - if not already existing:
    for (int i = 0; i < selecteddisjointregions.size(); i++)
    {
        int disjreg = selecteddisjointregions[i];
        
        // Get the element type number in the current disjoint region:
        int elementtypenumber = mydisjointregions->getelementtypenumber(disjreg);
        int elementdimension = mydisjointregions->getelementdimension(disjreg);

        // Get the number of form functions:
        std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, fieldtoadd->gettypename());
        int numberofformfunctions = myformfunction->count(fieldtoadd->getinterpolationorder(disjreg), elementdimension, 0);
        
        // Only treat the form functions not yet in the dof structure.
        if (rangebegin[fieldindex][disjreg].size() < numberofformfunctions)
        {
            int numffdefinedbeforeresize = rangebegin[fieldindex][disjreg].size();
            int currentnumberofdofs = mydisjointregions->countelements(disjreg);

            rangebegin[fieldindex][disjreg].resize(numberofformfunctions);
            rangeend[fieldindex][disjreg].resize(numberofformfunctions);
                
            for (int ff = numffdefinedbeforeresize; ff < numberofformfunctions; ff++)
            {
                rangebegin[fieldindex][disjreg][ff] = numberofdofs;
                rangeend[fieldindex][disjreg][ff] = numberofdofs + currentnumberofdofs - 1;
                numberofdofs += currentnumberofdofs;
            }
        }
    }
}

void dofmanager::addtostructure(shared_ptr<rawfield> fieldtoadd, int physicalregionnumber)
{
    // Get all disjoint regions in the physical region with (-1):
	std::vector<int> disjregs = ((universe::mymesh->getphysicalregions())->get(physicalregionnumber))->getdisjointregions(-1);
    
    addtostructure(fieldtoadd, disjregs);
}

void dofmanager::selectfield(shared_ptr<rawfield> selectedfield)
{
    selectedfieldnumber = -1;
    for (int i = 0; i < myfields.size(); i++)
    {
        if (myfields[i].get() == selectedfield.get())
        {
            selectedfieldnumber = i;
            break;
        }
    }
    if (selectedfieldnumber == -1)
    {
        std::cout << "Error in object 'dofmanager': selected field is not defined in the dof manager" << std::endl;
        abort();
    }
}

std::vector<int> dofmanager::getdisjointregionsofselectedfield(void)
{
    int totalnumdisjreg = universe::mymesh->getdisjointregions()->count();
    
    std::vector<int> output(totalnumdisjreg);
    int index = 0;
    for (int i = 0; i < totalnumdisjreg; i++)
    {
        if (rangebegin[selectedfieldnumber][i].size() > 0)
        {
            output[index] = i;
            index++;
        }
    }
    output.resize(index);
    
    return output;
}

int dofmanager::getrangebegin(int disjreg, int formfunc) { return rangebegin[selectedfieldnumber][disjreg][formfunc]; }
int dofmanager::getrangeend(int disjreg, int formfunc) { return rangeend[selectedfieldnumber][disjreg][formfunc]; }

int dofmanager::countconstraineddofs(void)
{
    int numconstraineddofs = 0;
    
    for (int fieldindex = 0; fieldindex < rangebegin.size(); fieldindex++)
    {
        for (int disjreg = 0; disjreg < rangebegin[fieldindex].size(); disjreg++)
        {
            // If the field is constrained on the disjoint region and there is at least one form function:
            if (rangebegin[fieldindex][disjreg].size() > 0 && myfields[fieldindex]->isconstrained(disjreg))
                numconstraineddofs += rangebegin[fieldindex][disjreg].size() * (rangeend[fieldindex][disjreg][0] - rangebegin[fieldindex][disjreg][0] + 1);
        }
    }
    return numconstraineddofs;
}

intdensematrix dofmanager::getconstrainedindexes(void)
{
    intdensematrix output(1,countconstraineddofs());
    int* myval = output.getvalues();
    
    int currentindex = 0;
    for (int fieldindex = 0; fieldindex < rangebegin.size(); fieldindex++)
    {
        for (int disjreg = 0; disjreg < rangebegin[fieldindex].size(); disjreg++)
        {
            // If the field is constrained on the disjoint region.
            if (myfields[fieldindex]->isconstrained(disjreg))
            {
                for (int ff = 0; ff < rangebegin[fieldindex][disjreg].size(); ff++)
                {
                    int numdofshere = rangeend[fieldindex][disjreg][0] - rangebegin[fieldindex][disjreg][0] + 1;
                    for (int i = 0; i < numdofshere; i++)
                    {
                        myval[currentindex] = rangebegin[fieldindex][disjreg][ff] + i;
                        currentindex++;
                    }
                }
            }
        }
    }
    return output;
}

void dofmanager::print(void)
{
    std::cout << "Showing the content of the dof manager (" << numberofdofs << " dofs in total):" << std::endl << std::endl;
    
    for (int i = 0; i < myfields.size(); i++)
    {
        for (int disjreg = 0; disjreg < rangebegin[i].size(); disjreg++)
        {
            // Only print not-empty entries:
            if (rangebegin[i][disjreg].size() != 0)
            {
                std::cout << "Field ";
                myfields[i]->print();
                std::cout << " @ disj. reg. " << disjreg << " is in range " << rangebegin[i][disjreg][0] << "-" << rangeend[i][disjreg][rangeend[i][disjreg].size()-1] << std::endl;
            }
        }
    }
    std::cout << std::endl;
}


intdensematrix dofmanager::getadresses(shared_ptr<rawfield> inputfield, int fieldinterpolationorder, int elementtypenumber, std::vector<int> &elementlist, int fieldphysreg, bool useminusonetag)
{
    elements* myelements = universe::mymesh->getelements();
    disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();

    // Create a vector whose ith index is true if the field is defined 
    // on the disjoint region number i and false otherwise.
    std::vector<int> fielddisjregs = ((universe::mymesh->getphysicalregions())->get(fieldphysreg))->getdisjointregions(-1); // Get all disj regs with (-1)
    std::vector<bool> isfielddefinedondisjointregion(mydisjointregions->count(),false);
    for (int i = 0; i < fielddisjregs.size(); i++)
        isfielddefinedondisjointregion[fielddisjregs[i]] = true;
                    
    element myelement(elementtypenumber);
    
    // Create a form function iterator to iterate through all form functions.
    hierarchicalformfunctioniterator myiterator(inputfield->gettypename(), elementtypenumber, fieldinterpolationorder);

    int numrows = myiterator.count();
    int numcols = elementlist.size();
    
    intdensematrix output(numrows, numcols);
    int* adresses = output.getvalues();

    selectfield(inputfield);
    
    for (int ff = 0; ff < myiterator.count(); ff++)
    {
        int associatedelementtype = myiterator.getassociatedelementtype();
        int formfunctionindex = myiterator.getformfunctionindexinnodeedgefacevolume();
        int num = myiterator.getnodeedgefacevolumeindex();
        // For quad subelements in prisms and pyramids:
        if ((elementtypenumber == 6 || elementtypenumber == 7) && associatedelementtype == 3)
            num -= myelement.counttriangularfaces();
        
        for (int i = 0; i < elementlist.size(); i++)
        {
            int elem = elementlist[i];
            // Get the subelement to which the current form function is associated:
            int currentsubelem = myelements->getsubelement(associatedelementtype, elementtypenumber, elem, num);
            // Also get its disjoint region number:
            int currentdisjointregion = myelements->getdisjointregion(associatedelementtype, currentsubelem);
            // If not in a disjoint region on which the field is defined set -2 adress.
            if (isfielddefinedondisjointregion[currentdisjointregion])
            {
                // Use it to get the subelem index in the disjoint region:
                currentsubelem -= mydisjointregions->getrangebegin(currentdisjointregion);

                adresses[ff*numcols+i] = rangebegin[selectedfieldnumber][currentdisjointregion][formfunctionindex] + currentsubelem;     

                if (useminusonetag && inputfield->isconstrained(currentdisjointregion))
                    adresses[ff*numcols+i] = -1;
            }
            else
                adresses[ff*numcols+i] = -2;
        }
        myiterator.next();
    }

    return output;
}
    



