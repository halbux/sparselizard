#include "oneconstant.h"
#include "element.h"

using namespace std;


oneconstant::oneconstant(int td, int et)
{
    targetdim = td;
    elementtypenumber = et;
    element myelem(elementtypenumber);
    elementdimension = myelem.getelementdimension();
}

int oneconstant::count(int order)
{
    if (order <= 0 || targetdim != elementdimension)
        return 0;
    else
        return 1;
}

int oneconstant::count(int order, int dim, int num)
{
    if (targetdim == elementdimension && dim == elementdimension)
        return count(order);
    else
        return 0;
}



hierarchicalformfunctioncontainer oneconstant::evalat(int maxorder) 
{    
    std::string type = "one"+std::to_string(targetdim);

    hierarchicalformfunctioncontainer val(type, elementtypenumber);

    polynomial formfunc;
    formfunc.set({{{1.0}}});
    val.set(1,elementdimension,0,0,0,0,formfunc);
    
    return val;
}

