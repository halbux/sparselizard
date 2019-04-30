#include "vectorfieldselect.h"


void vectorfieldselect::setdata(int physreg, field myfield, std::string op)
{
    (myfield.getpointer())->transferdata(physreg, vectorfieldselect(myrawvec, myrawfield), op);
}
