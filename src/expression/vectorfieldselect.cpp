#include "vectorfieldselect.h"


void vectorfieldselect::setdata(int physreg, field myfield, std::string op)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    
    if (op != "set" && op != "add")
    {
        logs log;
        log.msg() << "Error in 'vectorfieldselect' object: operation " << op << " is unknown in .setdata (use 'set' or 'add')" << std::endl;
        log.error();
    }
    
    (myfield.getpointer())->transferdata(physreg, vectorfieldselect(myrawvec, myrawfield), op);
}
