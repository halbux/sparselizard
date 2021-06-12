#include "vectorfieldselect.h"


void vectorfieldselect::setdata(int physreg, field myfield, std::string op)
{
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    
    if (op != "set" && op != "add")
    {
        std::cout << "Error in 'vectorfieldselect' object: operation " << op << " is unknown in .setdata (use 'set' or 'add')" << std::endl;
        abort();
    }
    
    (myfield.getpointer())->transferdata(physreg, vectorfieldselect(myrawvec, myrawfield), op);
}
