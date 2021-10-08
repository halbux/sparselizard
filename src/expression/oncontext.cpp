#include "oncontext.h"


oncontext::oncontext(int physreg, expression* coordshift, bool errorifnotfound)    
{
    myisdefined = true;

    myphysreg = physreg;
    myerrorifnotfound = errorifnotfound;
    if (coordshift == NULL)
        mycoordshift = {};
    else
        mycoordshift = {*coordshift};
}     

bool oncontext::isdefined(void)
{
    return myisdefined;
}  

int oncontext::getphysicalregion(void) 
{
    return myphysreg;
}

bool oncontext::isshifted(void)
{
    return (mycoordshift.size() != 0);
}

expression* oncontext::getshift(void)
{
    return &(mycoordshift[0]);
}

bool oncontext::iserrorifnotfound(void)
{
    return myerrorifnotfound;
}

bool oncontext::isequal(oncontext* tocompare)
{
    bool isitequal = true;
    
    if (myisdefined != tocompare->myisdefined)
        isitequal = false;
    
    if (myphysreg != tocompare->myphysreg)
        isitequal = false;
        
    if (myerrorifnotfound != tocompare->myerrorifnotfound)
        isitequal = false;
    
    // Compare the coordinate shift expressions:
    if (mycoordshift.size() != tocompare->mycoordshift.size())
        isitequal = false;
    else
    {
        if (mycoordshift.size() != 0)
        {
            expression cs1 = mycoordshift[0];
            expression cs2 = tocompare->mycoordshift[0];
        
            // Compare the two expressions:
            if (cs1.countrows() != cs2.countrows() || cs1.countcolumns() != cs2.countcolumns())
                isitequal = false;
            else
            {
                for (int i = 0; i < cs1.countrows(); i++)
                {
                    for (int j = 0; j < cs1.countcolumns(); j++)
                    {
                        if (cs1.getoperationinarray(i,j).get() != cs2.getoperationinarray(i,j).get())
                            isitequal = false;
                    }
                }
            }
        }
    }
    
    return isitequal;
}

