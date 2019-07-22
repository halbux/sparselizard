#include "oneconstant.h"

using namespace std;


int oneconstant::count(int order)
{
    if (order <= 0)
        return 0;
    
    element myelement(myelementtypenumber);
    int elemdim = myelement.getelementdimension();
    
	int problemdimension = universe::mymesh->getmeshdimension();
	
    if (elemdim == problemdimension)
    	return 1;
	else
		return 0;
}

int oneconstant::count(int order, int dim, int num)
{
    // The 'num' input argument is not required here since all nodes, 
    // edges and faces have the same number of form functions. It is
    // however required for prisms and pyramids.
    
    if (order <= 0)
        return 0;
    
    element myelement(myelementtypenumber);
    int elemdim = myelement.getelementdimension();
    
	int problemdimension = universe::mymesh->getmeshdimension();
    
    if (elemdim == problemdimension && dim == problemdimension)
    	return 1;
	else
		return 0;
}



hierarchicalformfunctioncontainer oneconstant::evalat(int maxorder) 
{    
    element myelement(myelementtypenumber);
    int elemdim = myelement.getelementdimension();
    
	int problemdimension = universe::mymesh->getmeshdimension();

    hierarchicalformfunctioncontainer val("one", myelementtypenumber);

	if (elemdim == problemdimension)
	{
		// Loop on all possible orientations. The value is the same.
		for (int orientation = 0; orientation < orientation::countorientations(myelementtypenumber); orientation++)
		{
			polynomial formfunc;
			formfunc.set({{{1.0}}});
			val.set(1,problemdimension,0,orientation,0,0,formfunc);
		}
	}
    
    return val;
}
