#include "hierarchicalformfunctioniterator.h"


hierarchicalformfunctioniterator::hierarchicalformfunctioniterator(std::string formfunctiontypename, int elementtypenumber, int order)
{
    myformfunctiontypename = formfunctiontypename;
    myelementtypenumber = elementtypenumber;
    myorder = order;
    
    // Set the current form function to the first one.
    next();
}

int hierarchicalformfunctioniterator::count(void)
{
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(myelementtypenumber, myformfunctiontypename);
    return myformfunction->count(myorder);
}

void hierarchicalformfunctioniterator::next()
{
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(myelementtypenumber, myformfunctiontypename);

    // Check if there is any other form function on the current node/edge/face/volume.
    // In case 'next' is called for the first time currentformfunctionindexinnodeedgefacevolume is at -1.
    if (currentformfunctionindexinnodeedgefacevolume == -1 || currentformfunctionindexinnodeedgefacevolume < myformfunction->count(myorder, currentdimension, currentnodeedgefacevolumeindex) - 1)
        currentformfunctionindexinnodeedgefacevolume++;
    else
    {
        // Otherwise go to the next node/edge/face/volume, if any:
        element myelement(myelementtypenumber);
        if (currentnodeedgefacevolumeindex < myelement.countdim(currentdimension) - 1)
        {
            currentnodeedgefacevolumeindex++;
            currentformfunctionindexinnodeedgefacevolume = 0;
        }
        else
        {
            // Otherwise go to the next dimension, if any:
            if (currentdimension < 3)
            {
                currentdimension++;
                currentnodeedgefacevolumeindex = 0;
                currentformfunctionindexinnodeedgefacevolume = 0;
            }
            // In case there are no more form functions:
            else
                return;
        }
    }
    
    // Make sure the form function is defined:
    if (currentformfunctionindexinnodeedgefacevolume >= myformfunction->count(myorder, currentdimension, currentnodeedgefacevolumeindex))
        next();
    else
        overallformfunctionindex++;
        
}

int hierarchicalformfunctioniterator::getformfunctionindexincurrentorderinnodeedgefacevolume(void)
{
    int currentorder = getformfunctionorder();
    
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(myelementtypenumber, myformfunctiontypename);

    return currentformfunctionindexinnodeedgefacevolume - myformfunction->count(currentorder - 1, currentdimension, currentnodeedgefacevolumeindex);
}

int hierarchicalformfunctioniterator::getformfunctionorder(void)
{
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(myelementtypenumber, myformfunctiontypename);
    
    for (int i = 0; i <= myorder; i++)
    {
        if (currentformfunctionindexinnodeedgefacevolume < myformfunction->count(i, currentdimension, currentnodeedgefacevolumeindex) )
            return i;
    }
    
    abort(); // fix return warning
}

int hierarchicalformfunctioniterator::getassociatedelementtype(void)
{
    element myelement(myelementtypenumber);

    switch (currentdimension)
    {
        case 0:
            // It is a point.
            return 0;
        case 1:
            // It is a line:
            return 1;
        case 2:
            // It is a triangle or a quadrangle:
            if (myelement.istriangularface(currentnodeedgefacevolumeindex))
                return 2;
            return 3;
        case 3:
            // It is myelementtypenumber itself!
            return myelementtypenumber;
    }
    
    abort(); // fix return warning
}

void hierarchicalformfunctioniterator::print(void)
{
    std::cout << "This is form function number " << overallformfunctionindex << ", associated to ";
    switch (currentdimension)
    {
        case 0:
            std::cout << "vertex";
            break;
        case 1:
            std::cout << "edge";
            break;
        case 2:
            std::cout << "face";
            break;
        case 3:
            std::cout << "volume";
            break;
    }
    std::cout << " number " << currentnodeedgefacevolumeindex << " and number " << currentformfunctionindexinnodeedgefacevolume << " of its kind (order is " << getformfunctionorder() << ")" << std::endl;
}


